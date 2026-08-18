#!/bin/bash -l
# Merge the per-array resource_usage TSVs into one table, enriched with Slurm accounting.
#
# The in-job tracker (src/common/resource_tracking.sh) records what /usr/bin/time measured:
# exact cumulative CPU, exact wall clock, and the peak RSS of the largest single process.
# That last figure undercounts methods that fork parallel workers, because ru_maxrss is a
# per-process high-water mark rather than the sum of concurrently resident processes.
#
# Slurm samples the whole process tree instead (JobAcctGatherType=jobacct_gather/linux, every
# JobAcctGatherFrequency=30s here), so its MaxRSS covers concurrent workers but can miss a
# spike shorter than the sampling interval. The two measures bracket real peak memory, which
# is why this script keeps both rather than picking one:
#
#   max_rss_kb        LOWER bound: ru_maxrss is the largest single process, so concurrent
#                     workers are missed entirely.
#   sacct_max_rss_kb  UPPER bound: jobacct_gather/linux sums RSS over the process tree, and
#                     forked workers sharing copy-on-write pages have those pages counted
#                     once per process. Also sampled every 30s, so short spikes are missed.
#
# So the true peak lies between them, and neither is "the" answer. Measured here: for
# single-process methods (Pando, CellOracle) the two agree within 7%, while for FigR -- which
# forks `future` workers -- sacct reads 4.2x higher than time -v.
#
# pss_peak_kb supersedes both where present. The tracker now samples the process tree's
# proportional set size directly, which divides each shared page by the number of processes
# mapping it and so counts genuinely shared memory exactly once. Validated against a synthetic
# workload holding 4 GB across 1 parent + 6 forked readers: PSS reported 4.01 GB (correct)
# while the summed-RSS reading was 28.05 GB, inflated exactly 7x. Prefer MED_PSS_GB; the LO/HI
# bracket remains for runs recorded before the sampler existed.
#
# NOTE ON THIS CLUSTER: cgroup.conf sets ConstrainCores=yes but does NOT set
# ConstrainRAMSpace, so --mem is a scheduling reservation only and is not enforced. Methods
# can and do exceed it (FigR and SCENIC+ rows above 72G on a 72G allocation, none OOM-killed).
# Equal-resource comparisons therefore hold for cores but not for memory.
#
# Usage:
#   src/common/collect_resource_usage.sh [-o OUT.tsv] [record.tsv ...]
#
# With no file arguments it collects every LOGS/run_tools/run_tools_*/resource_usage_*.tsv.
# Writes to stdout unless -o is given.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUT=""

while getopts ":o:h" opt; do
    case "$opt" in
        o) OUT="$OPTARG" ;;
        h) sed -n '2,35p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# ===== INPUT RECORDS =====
records=("$@")
if [ "${#records[@]}" -eq 0 ]; then
    # nullglob so a project with no runs yet produces a clean error instead of a literal glob.
    shopt -s nullglob
    records=("${PROJECT_DIR}"/LOGS/run_tools/run_tools_*/resource_usage_*.tsv)
    shopt -u nullglob
fi

if [ "${#records[@]}" -eq 0 ]; then
    echo "ERROR: no resource_usage_*.tsv found under ${PROJECT_DIR}/LOGS/run_tools/" >&2
    exit 1
fi

echo "Collecting from ${#records[@]} record file(s)" >&2

# ===== SLURM ACCOUNTING =====
# One sacct call for every job id mentioned, rather than one per row.
job_ids=$(awk -F'\t' 'FNR > 1 && $16 != "NA" { print $16 }' "${records[@]}" | sort -u | paste -sd,)

sacct_raw=$(mktemp)
trap 'rm -f "${sacct_raw}"' EXIT

if [ -n "${job_ids}" ]; then
    # Two rows matter per task: the job-level row carries State and ReqMem, while the .batch
    # step is the only one with MaxRSS (it accounts for the batch script and all its
    # children; .extern covers just Slurm's own bookkeeping). Both are parsed below.
    # sacct -P delimits with '|', so translate to tabs to match the record files and let one
    # awk read both with a single FS.
    sacct -j "${job_ids}" \
        --format=JobID,State,Elapsed,TotalCPU,MaxRSS,MaxVMSize,ReqMem,NCPUS \
        -P -n 2>/dev/null | tr '|' '\t' > "${sacct_raw}" || true
fi

if [ ! -s "${sacct_raw}" ]; then
    echo "WARNING: sacct returned nothing; sacct_* columns will be NA" >&2
fi

# ===== JOIN =====
emit() {
    awk -F'\t' -v OFS='\t' '
        # "4998964K", "1.20G", "523052K" -> kilobytes. Slurm suffixes K/M/G/T (and sometimes
        # none, which it means as kilobytes).
        function to_kb(v,   n, u) {
            if (v == "" || v == "NA") return "NA"
            u = substr(v, length(v), 1)
            n = v + 0
            if (u == "K" || u == "k") return sprintf("%.0f", n)
            if (u == "M" || u == "m") return sprintf("%.0f", n * 1024)
            if (u == "G" || u == "g") return sprintf("%.0f", n * 1048576)
            if (u == "T" || u == "t") return sprintf("%.0f", n * 1073741824)
            return sprintf("%.0f", n)
        }
        # "1-02:03:04", "05:07:12", "13.124", "00:00:00" -> seconds.
        function to_s(v,   days, rest, p, n, i, s) {
            if (v == "" || v == "NA") return "NA"
            days = 0; rest = v
            if (index(v, "-") > 0) { split(v, p, "-"); days = p[1] + 0; rest = p[2] }
            n = split(rest, p, ":")
            s = 0
            for (i = 1; i <= n; i++) s = s * 60 + p[i] + 0
            return sprintf("%.2f", days * 86400 + s)
        }

        # --- pass 1: sacct rows ---
        FNR == NR {
            id = $1
            if (id ~ /\.extern$/) next
            if (id ~ /\.batch$/) {
                # The step row: the only place MaxRSS/MaxVMSize are populated.
                sub(/\.batch$/, "", id)
                elapsed[id] = to_s($3)
                cpu[id]     = to_s($4)
                rss[id]     = to_kb($5)
                vmem[id]    = to_kb($6)
            } else if (id !~ /\./) {
                # The job row: authoritative State, and the only place ReqMem appears.
                state[id]  = $2
                reqmem[id] = $7
                ncpus[id]  = $8
                # Kept as a fallback for jobs whose .batch row has been purged.
                if (!(id in elapsed)) elapsed[id] = to_s($3)
                if (!(id in cpu))     cpu[id]     = to_s($4)
            }
            next
        }

        # --- pass 2: tracker records ---
        # Header: emit once, extended with the sacct columns.
        FNR == 1 {
            if (!header_done) {
                print $0, "sacct_state", "sacct_elapsed_s", "sacct_total_cpu_s",
                          "sacct_max_rss_kb", "sacct_max_rss_gb", "sacct_max_vmsize_kb",
                          "sacct_req_mem", "sacct_ncpus"
                header_done = 1
            }
            next
        }
        {
            # slurm_job_id is column 16, array_task_id 17 (see RT_COLUMNS).
            key = ($17 == "NA" || $17 == "") ? $16 : $16 "_" $17
            r = (key in rss) ? rss[key] : "NA"
            print $0,
                  (key in state)   ? state[key]   : "NA",
                  (key in elapsed) ? elapsed[key] : "NA",
                  (key in cpu)     ? cpu[key]     : "NA",
                  r,
                  (r == "NA") ? "NA" : sprintf("%.3f", r / 1048576),
                  (key in vmem)    ? vmem[key]    : "NA",
                  (key in reqmem)  ? reqmem[key]  : "NA",
                  (key in ncpus)   ? ncpus[key]   : "NA"
        }
    ' "${sacct_raw}" "${records[@]}"
}

if [ -n "${OUT}" ]; then
    mkdir -p "$(dirname "${OUT}")"
    emit > "${OUT}"
    echo "Wrote ${OUT} ($(( $(wc -l < "${OUT}") - 1 )) rows)" >&2

    # ===== PER-METHOD SUMMARY =====
    # Median is more informative than mean here: a handful of runs die in seconds on a bad
    # node or a stale R library, and a mean over those is meaningless.
    echo >&2
    awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        $c["step"] != "TOTAL" || $c["exit_status"] != "0" { next }
        {
            m = $c["method"]
            n[m]++
            wall[m, n[m]] = $c["wall_s"] + 0
            cpu[m, n[m]]  = $c["cpu_s"] + 0
            # Both bounds are carried separately; see the header comment on why collapsing
            # them to one number misranks forking methods against single-process ones.
            rss_lo[m, n[m]] = ($c["max_rss_gb"] == "NA") ? 0 : $c["max_rss_gb"] + 0
            rss_hi[m, n[m]] = ($c["sacct_max_rss_gb"] == "NA") ? 0 : $c["sacct_max_rss_gb"] + 0
            # PSS is the measured truth where the sampler ran; kept in its own tally so runs
            # predating it do not drag a median toward zero.
            if ("pss_peak_gb" in c && $c["pss_peak_gb"] != "NA") {
                P[m, ++np[m]] = $c["pss_peak_gb"] + 0
            }
        }
        function median(m, arr, cnt,   v, i, j, t) {
            for (i = 1; i <= cnt; i++) v[i] = arr[m, i]
            for (i = 1; i < cnt; i++) for (j = i + 1; j <= cnt; j++)
                if (v[j] < v[i]) { t = v[i]; v[i] = v[j]; v[j] = t }
            return (cnt % 2) ? v[(cnt + 1) / 2] : (v[cnt / 2] + v[cnt / 2 + 1]) / 2
        }
        END {
            printf "%-12s %6s  %12s  %12s  %12s  %12s  %12s\n", "METHOD", "RUNS", "MED_WALL_H", "MED_CPU_H", "MED_PSS_GB", "MED_RSS_LO", "MED_RSS_HI"
            for (m in n) {
                pss = (np[m] > 0) ? sprintf("%12.2f", median(m, P, np[m])) : sprintf("%12s", "-")
                printf "%-12s %6d  %12.2f  %12.2f  %s  %12.2f  %12.2f\n", m, n[m],
                    median(m, wall, n[m]) / 3600, median(m, cpu, n[m]) / 3600,
                    pss, median(m, rss_lo, n[m]), median(m, rss_hi, n[m])
            }
        }
    ' "${OUT}" | { IFS= read -r hdr; printf '%s\n' "${hdr}"; sort -k1,1; } >&2
else
    emit
fi
