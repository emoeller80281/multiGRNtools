#!/bin/bash
# Resource tracking for GRN inference methods.
#
# The method scripts under src/<Method>/ already wrap their individual steps in
# `/usr/bin/time -v`, which tells us what each step cost but never what the method cost as
# a whole. This library adds that missing total: rt_run_method wraps one entire method
# invocation in a single `/usr/bin/time -v`, so the numbers cover every child process the
# method spawns.
#
# Why one outer wrapper rather than summing the per-step logs: GNU time reports
# getrusage(RUSAGE_CHILDREN), which the kernel accumulates over *all* waited-for
# descendants. So a single outer wrapper already gives cumulative CPU across every step and
# subprocess, plus true end-to-end wall time. Summing per-step logs would instead miss
# everything the method does between steps (conda activation, file staging, formatting).
#
# The one number the outer wrapper does NOT give cleanly is aggregate peak memory: ru_maxrss
# is the high-water mark of the single largest process, not the sum of concurrently resident
# ones. Methods that fork parallel workers (FigR/Pando via future, LINGER via num_cpu,
# SCENIC+ via snakemake) can therefore hold more than max_rss_kb reports.
#
# collect_resource_usage.sh joins Slurm's own MaxRSS onto these records as a second estimate.
# The two bracket the true peak rather than either being correct:
#   max_rss_kb        LOWER bound -- largest single process, concurrent workers missed
#   sacct_max_rss_kb  UPPER bound -- RSS summed over the process tree, so pages shared by
#                     forked workers are counted once per process; also sampled every 30s
# Do not collapse them to one number: measured on this project the two agree within 7% for
# single-process methods but differ 4.2x for FigR, so picking either alone reorders the
# methods relative to each other.
#
# Usage (see run_tools_stability.sh / run_tools.sh):
#   source "${PROJECT_DIR}/src/common/resource_tracking.sh"
#   rt_init                                  # picks up RT_RECORD_FILE, sets defaults
#   rt_run_method "$METHOD" "$script" "$stdout_log" "$stderr_log"
#
# Records are appended to $RT_RECORD_FILE as one tab-separated line per method run. Each
# line is written with a single printf so concurrent array tasks appending to the same file
# cannot interleave (an O_APPEND write below PIPE_BUF is atomic) -- the same approach
# run_tools_stability.sh already uses for its failures file.

# Column order of $RT_RECORD_FILE. Kept in one place so the collector can assert on it.
RT_COLUMNS=(
    method cell_type sample_name subsample step
    wall_s user_s sys_s cpu_s cpu_pct
    max_rss_kb max_rss_gb fs_inputs fs_outputs
    exit_status slurm_job_id array_task_id node ncpus started_at
    pss_peak_kb pss_peak_gb rss_sum_peak_kb pss_samples
)

# How often the PSS sampler walks the process tree, in seconds. Reading smaps_rollup costs
# ~0.4 ms per process even at 8 GB, so 10s is far cheaper than it looks; raise it only if a
# method spawns thousands of short-lived processes.
RT_PSS_INTERVAL="${RT_PSS_INTERVAL:-10}"

# Absolute path to GNU time. Only the compute nodes have it (the login nodes do not ship the
# `time` RPM), so this is resolved at call time and the tracker degrades to wall-clock-only
# rather than failing the run.
RT_TIME_BIN="${RT_TIME_BIN:-/usr/bin/time}"

rt_init() {
    # Where the parsed one-line-per-run records go. Callers normally set this to sit next to
    # the job's other logs; the fallback keeps the library usable when sourced ad hoc.
    RT_RECORD_FILE="${RT_RECORD_FILE:-${LOG_DIR:-.}/resource_usage.tsv}"
    # Where the raw `/usr/bin/time -v` reports go, one file per method run, kept for the
    # fields this library does not parse (context switches, page faults, VM size).
    RT_RAW_DIR="${RT_RAW_DIR:-$(dirname "${RT_RECORD_FILE}")/resource_raw}"

    mkdir -p "$(dirname "${RT_RECORD_FILE}")" "${RT_RAW_DIR}"

    # Header once per file. Written under a noclobber redirect so that when 100 array tasks
    # race here, exactly one wins and the rest silently skip instead of appending duplicate
    # or truncating the file.
    if [ ! -s "${RT_RECORD_FILE}" ]; then
        local header
        header=$(printf '%s\t' "${RT_COLUMNS[@]}")
        ( set -o noclobber; printf '%s\n' "${header%$'\t'}" > "${RT_RECORD_FILE}" ) 2>/dev/null || true
    fi
}

# rt_parse_time_v <raw_time_v_file>
# Turns a GNU `time -v` report into the tab-separated middle of a record:
#   wall_s user_s sys_s cpu_s cpu_pct max_rss_kb max_rss_gb fs_inputs fs_outputs exit_status
# Emits "NA" for anything the report did not contain, so the column count is always fixed.
rt_parse_time_v() {
    local f="$1"
    [ -s "$f" ] || { printf 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA'; return 0; }

    awk -F': ' '
        # Elapsed comes formatted as h:mm:ss or m:ss.ss, so convert to plain seconds.
        function to_seconds(t,   n, p, s) {
            n = split(t, p, ":")
            s = 0
            for (i = 1; i <= n; i++) s = s * 60 + p[i] + 0
            return s
        }
        /User time \(seconds\)/                 { user = $2 + 0 }
        /System time \(seconds\)/               { sys  = $2 + 0 }
        /Percent of CPU this job got/           { pct  = $2; sub(/%/, "", pct); pct = pct + 0 }
        /Elapsed \(wall clock\) time/           { wall = to_seconds($NF) }
        /Maximum resident set size \(kbytes\)/  { rss  = $2 + 0 }
        /File system inputs/                    { fsin = $2 + 0 }
        /File system outputs/                   { fsout= $2 + 0 }
        # A killed child is reported as "Command terminated by signal N" on the FIRST line and
        # still gets a trailing "Exit status: 0", so the signal has to win -- otherwise a
        # SIGKILLed method is recorded as a clean exit. Surfaced as 128+N the way a shell does.
        /Exit status/                           { if (!signalled) code = $2 + 0 }
        /Command terminated by signal/          { split($0, s, "signal "); code = 128 + (s[2] + 0); signalled = 1 }
        END {
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                (wall == "" ? "NA" : sprintf("%.2f", wall)),
                (user == "" ? "NA" : sprintf("%.2f", user)),
                (sys  == "" ? "NA" : sprintf("%.2f", sys)),
                (user == "" || sys == "" ? "NA" : sprintf("%.2f", user + sys)),
                (pct  == "" ? "NA" : pct),
                (rss  == "" ? "NA" : rss),
                (rss  == "" ? "NA" : sprintf("%.3f", rss / 1048576)),
                (fsin == "" ? "NA" : fsin),
                (fsout== "" ? "NA" : fsout),
                (code == "" ? "NA" : code)
        }
    ' "$f"
}

# rt_pss_sample <root_pid>
# Prints "<pss_kb> <rss_sum_kb>" for the whole process tree rooted at <root_pid>.
#
# Pss (proportional set size) divides each shared page by the number of processes mapping it,
# so summing Pss across a tree counts genuinely shared memory exactly once. That is the number
# neither existing measure gives: ru_maxrss sees only the largest single process, and Slurm
# sums Rss so copy-on-write pages are counted once per forked worker.
#
# rss_sum is summed alongside from the same samples purely as a diagnostic: it reproduces what
# Slurm's accounting does, so rss_sum/pss is a direct read-out of how much a method's memory is
# shared rather than real.
rt_pss_sample() {
    ps -eo pid=,ppid= 2>/dev/null | awk -v root="$1" '''
        { ppid[$1] = $2; pid[++n] = $1 }
        END {
            # Walk the parent map to a fixpoint to collect every descendant of root.
            desc[root] = 1
            do {
                changed = 0
                for (i = 1; i <= n; i++) {
                    p = pid[i]
                    if (!(p in desc) && (ppid[p] in desc)) { desc[p] = 1; changed = 1 }
                }
            } while (changed)

            pss = 0; rss = 0
            for (p in desc) {
                f = "/proc/" p "/smaps_rollup"
                # A process that exits mid-walk simply yields no lines; that is fine, the next
                # sample picks up whatever is alive then.
                while ((getline line < f) > 0) {
                    if (line ~ /^Pss:/)      { split(line, a, /[ \t]+/); pss += a[2] }
                    else if (line ~ /^Rss:/) { split(line, a, /[ \t]+/); rss += a[2] }
                }
                close(f)
            }
            printf "%d %d\n", pss, rss
        }'''
}

# rt_pss_watch <root_pid> <outfile> <interval_s>
# Polls until <root_pid> exits, keeping the running peak in <outfile> as
# "<peak_pss_kb> <peak_rss_sum_kb> <n_samples>". Meant to be backgrounded.
rt_pss_watch() {
    local root="$1" out="$2" interval="$3"
    local peak_pss=0 peak_rss=0 samples=0 cur_pss cur_rss line

    while kill -0 "$root" 2>/dev/null; do
        line=$(rt_pss_sample "$root")
        cur_pss=${line%% *}
        cur_rss=${line##* }
        samples=$((samples + 1))
        if [ "${cur_pss:-0}" -gt "$peak_pss" ]; then peak_pss=$cur_pss; fi
        if [ "${cur_rss:-0}" -gt "$peak_rss" ]; then peak_rss=$cur_rss; fi
        printf '%d %d %d\n' "$peak_pss" "$peak_rss" "$samples" > "$out"
        sleep "$interval"
    done
}

# rt_read_pss <pssfile>
# Renders the sampler's peaks as the four trailing record fields, or NAs if it never ran
# (no /proc/smaps_rollup, or a method that finished inside the first interval).
rt_read_pss() {
    local f="$1"
    if [ -s "$f" ]; then
        awk '{ printf "%d\t%.3f\t%d\t%d", $1, $1 / 1048576, $2, $3 }' "$f"
    else
        printf 'NA\tNA\tNA\tNA'
    fi
}

# rt_slug
# Identifies this particular run for use in filenames. RT_RAW_DIR is shared by every task in
# an array, and two tasks routinely run the same method on different subsamples, so the raw
# report name has to carry the sample identity or they overwrite each other.
rt_slug() {
    local slug="${CELL_TYPE:-na}_${SAMPLE_NAME:-na}"
    [ -n "${SUBSAMPLE_NUM:-}" ] && slug="${slug}_sub${SUBSAMPLE_NUM}"
    [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && slug="${slug}_task${SLURM_ARRAY_TASK_ID}"
    # Collapse anything that would be awkward in a filename.
    printf '%s' "${slug//[^A-Za-z0-9._-]/_}"
}

# rt_record <method> <step> <parsed_fields_tsv> <started_at>
# Appends one record, prefixing the run's identity and suffixing the Slurm context. Also
# leaves the line in RT_LAST_RECORD so the caller can report its own numbers without having
# to re-read a file that every other array task is appending to.
rt_record() {
    local method="$1" step="$2" parsed="$3" started_at="$4" pss="${5:-$(printf 'NA\tNA\tNA\tNA')}"

    RT_LAST_RECORD=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "${method}" \
        "${CELL_TYPE:-NA}" \
        "${SAMPLE_NAME:-NA}" \
        "${SUBSAMPLE_NUM:-NA}" \
        "${step}" \
        "${parsed}" \
        "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-NA}}" \
        "${SLURM_ARRAY_TASK_ID:-NA}" \
        "${SLURMD_NODENAME:-$(hostname -s)}" \
        "${SLURM_CPUS_PER_TASK:-NA}" \
        "${started_at}" \
        "${pss}")

    printf '%s\n' "${RT_LAST_RECORD}" >> "${RT_RECORD_FILE}"
}

# rt_run_method <method> <script> <stdout_log> <stderr_log>
# Runs one method end to end under `/usr/bin/time -v` and records its total cost.
# Returns the method's own exit status, so callers under `set -e` keep aborting exactly as
# they did before this tracking existed.
rt_run_method() {
    local method="$1" script="$2" stdout_log="$3" stderr_log="$4"
    local started_at raw rc=0

    started_at=$(date +%Y-%m-%dT%H:%M:%S)
    raw="${RT_RAW_DIR}/${method}_$(rt_slug).time"
    RT_LAST_RAW="${raw}"

    local pss_file="${RT_RAW_DIR}/${method}_$(rt_slug).pss"
    rm -f "${pss_file}"

    if [ -x "${RT_TIME_BIN}" ]; then
        # Backgrounded so the PSS sampler has a pid to follow; `wait` below still yields the
        # method's own exit status, so callers under set -e behave exactly as before.
        # -o keeps the report out of the method's own stderr, so stderr_log stays clean for
        # the method's messages and the report survives even when the method fails.
        "${RT_TIME_BIN}" -v -o "${raw}" \
            bash "${script}" > "${stdout_log}" 2> "${stderr_log}" &
        local method_pid=$!
        rt_pss_watch "${method_pid}" "${pss_file}" "${RT_PSS_INTERVAL}" &
        local watch_pid=$!

        wait "${method_pid}" || rc=$?
        # The watcher exits on its own once the method is gone, but kill it explicitly so a
        # method that dies during a sleep does not hold the job open for one more interval.
        kill "${watch_pid}" 2>/dev/null || true
        wait "${watch_pid}" 2>/dev/null || true

        rt_record "${method}" "TOTAL" "$(rt_parse_time_v "${raw}")" "${started_at}" \
            "$(rt_read_pss "${pss_file}")"
    else
        # No GNU time on this host: still record wall clock and exit status so the row exists
        # and the missing columns are explicitly NA rather than silently absent.
        local t0=${SECONDS}
        bash "${script}" > "${stdout_log}" 2> "${stderr_log}" &
        local method_pid=$!
        rt_pss_watch "${method_pid}" "${pss_file}" "${RT_PSS_INTERVAL}" &
        local watch_pid=$!
        wait "${method_pid}" || rc=$?
        kill "${watch_pid}" 2>/dev/null || true
        wait "${watch_pid}" 2>/dev/null || true
        local wall=$(( SECONDS - t0 ))
        printf 'WARNING: %s not found, recording wall clock only\n' "${RT_TIME_BIN}" >&2
        rt_record "${method}" "TOTAL" \
            "$(printf '%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%d' "${wall}" "${rc}")" \
            "${started_at}" "$(rt_read_pss "${pss_file}")"
    fi

    return ${rc}
}

# rt_summary_last
# Human-readable report of the run this shell just performed, so the totals show up in the
# Slurm .out next to everything else about the task. Reads RT_LAST_RECORD rather than
# RT_RECORD_FILE: the file is shared by the whole array, so its last line usually belongs to
# a different task.
rt_summary_last() {
    [ -n "${RT_LAST_RECORD:-}" ] || return 0

    # "${RT_COLUMNS[*]}" joins on the first character of IFS, so pin IFS to a space here.
    # Callers legitimately leave IFS set to something else (run_tools*.sh parse the pipe-
    # delimited EXPERIMENT_LIST), and without this the awk below sees one giant column name
    # and every lookup silently resolves to field 1.
    local IFS=' '

    printf '%s\n' "${RT_LAST_RECORD}" | awk -F'\t' -v cols="${RT_COLUMNS[*]}" '
        BEGIN { split(cols, name, " "); for (i in name) col[name[i]] = i }
        {
            printf "  method    : %s\n", $col["method"]
            printf "  run       : %s/%s/subsample_%s\n", $col["cell_type"], $col["sample_name"], $col["subsample"]
            printf "  wall time : %s s (%.2f h)\n", $col["wall_s"], $col["wall_s"] / 3600
            printf "  cpu time  : %s s (user %s + sys %s), %s%% of one core\n", $col["cpu_s"], $col["user_s"], $col["sys_s"], $col["cpu_pct"]
            printf "  peak rss  : %s kB (%s GB, largest single process)\n", $col["max_rss_kb"], $col["max_rss_gb"]
            if ($col["pss_peak_kb"] != "NA") {
                printf "  peak pss  : %s kB (%s GB, whole tree, shared pages counted once)\n", $col["pss_peak_kb"], $col["pss_peak_gb"]
                printf "  rss sum   : %s kB (same samples, shared pages counted per-process)\n", $col["rss_sum_peak_kb"]
                printf "  samples   : %s at %ss\n", $col["pss_samples"], ENVIRON["RT_PSS_INTERVAL"]
            }
            printf "  file i/o  : %s in / %s out\n", $col["fs_inputs"], $col["fs_outputs"]
            printf "  exit      : %s\n", $col["exit_status"]
        }
    '
}
