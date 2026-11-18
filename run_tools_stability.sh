#!/bin/bash -l
#SBATCH --job-name="multiGRNtools_stability"
#SBATCH --output="LOGS/run_tools_stability/run_tools_stability_%A/job_%a.out"
#SBATCH --error="LOGS/run_tools_stability/run_tools_stability_%A/job_%a.err"
#SBATCH --time=72:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 12
#SBATCH --mem=128G
#SBATCH --exclude=psh01com1ochm[01-03]
#SBATCH --array=0-299%30

set -euo pipefail

source /gpfs/Home/esm5360/miniconda3/etc/profile.d/conda.sh

# ===== SELECTED METHODS =====
METHOD_LIST=(
    "Pando"
    "FigR"
    "SCENIC_PLUS"
    "LINGER"
    "CellOracle"
)

# ===== SAMPLE CONFIGURATION =====
# tissue|sample_name|cell_type|species -- cell_type is the per-cell-type subdirectory under
# ${RAW_DATA_DIR}/<tissue>/<sample_name>/. Ordered largest cell count first so the longest
# runs are the first ones the array starts.
EXPERIMENT_LIST=(
    "mouse_liver|liver_sample|Hepatocytes|mouse"                        # 15,937 cells
    "GSE140203_shareseq|skin_late_anagen|Basal_keratinocytes|mouse"     # 10,447 cells
    "kidney|Ctrl_4weeks_2|Fibroblasts|mouse"                            #    415 cells
    "10x_E18_mouse_brain|E18_brain|Excitatory_neurons|mouse"            #  2,467 cells
    "mESC|E7.5_rep2|Neurectoderm|mouse"                                 #    998 cells
    "GSE246464_HSC|young_rep1|HSC|mouse"                                #  3,658 cells
)

# Contains the cell columns to use for each subsample
#{
#    "HSC": {
#        "dataset": "GSE246464_HSC",
#        "cell_type": "HSC",
#        "sample_name": "young_rep1",
#        "subsamples": {
#            "subsample_0": [
#                "ACCGGTTCAGAACCGA-1",
#                "ACTCAGTAGCTTAGTA-1",
#                "CCACAATAGGCAAGTA-1",
#                "CGCTACTTCCGTAAAC-1",
#                "GCTATCCTCAACAAGG-1",
#            ]
#        }
#    },
#	"Liver": {
#	"dataset": "mouse_liver",
#	"cell_type": "Hepatocytes",
#	"sample_name": "liver_sample",
#	"subsamples": {
#		"subsample_0": [
#			"TGTTAGCAGGCGCTAC-4",
#			"CATGCGGGTATTGTCC-4",
#			"GGTCAAGCATTGTCAG-4",
#		]
#	}
#}

STABILITY_JSON_FILE="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.SINGLE_CELL_GRN_INFERENCE.MOELLER/subsample_cell_dict.json"

NUM_EXPERIMENTS=${#EXPERIMENT_LIST[@]}
NUM_METHODS=${#METHOD_LIST[@]}
NUM_SUBSAMPLES=10
TOTAL_TASKS=$((NUM_EXPERIMENTS * NUM_METHODS * NUM_SUBSAMPLES))

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ "$TASK_ID" -ge "$TOTAL_TASKS" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${TASK_ID} exceeds TOTAL_TASKS=${TOTAL_TASKS}"
    exit 1
fi

# Experiment is the outermost index so the largest samples start first; within one
# experiment the methods vary fastest, then the subsample.
TASKS_PER_EXPERIMENT=$((NUM_METHODS * NUM_SUBSAMPLES))
EXPERIMENT_ID=$((TASK_ID / TASKS_PER_EXPERIMENT))
SUBSAMPLE_NUM=$(((TASK_ID % TASKS_PER_EXPERIMENT) / NUM_METHODS))
METHOD_ID=$((TASK_ID % NUM_METHODS))

EXPERIMENT_CONFIG="${EXPERIMENT_LIST[$EXPERIMENT_ID]}"
METHOD="${METHOD_LIST[$METHOD_ID]}"

IFS='|' read -r TISSUE SAMPLE_NAME CELL_TYPE SPECIES <<< "$EXPERIMENT_CONFIG"

# The method scripts name their formatted GRN <prefix>_${CELL_TYPE}_${SAMPLE_NAME}.tsv and key
# their working dirs on SAMPLE_NAME, so they are given this suffixed name. That gives each
# subsample its own GRN file and keeps concurrent subsamples of one sample apart.
RUN_SAMPLE_NAME="${SAMPLE_NAME}_sub${SUBSAMPLE_NUM}"

echo "TASK_ID=${TASK_ID}"
echo "EXPERIMENT_ID=${EXPERIMENT_ID}"
echo "METHOD_ID=${METHOD_ID}"
echo "METHOD=${METHOD}"
echo "TISSUE=${TISSUE}"
echo "SAMPLE_NAME=${SAMPLE_NAME}"
echo "CELL_TYPE=${CELL_TYPE}"
echo "SPECIES=${SPECIES}"
echo "SUBSAMPLE_NUM=${SUBSAMPLE_NUM}"


# ===== PATH CONFIGURATION =====
PROJECT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
RAW_DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/MUON_FILTERED_COUNT_DATASETS/mm10"
DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/data"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
REFERENCE_GENOME_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.SINGLE_CELL_GRN_INFERENCE.MOELLER/data/genome_data/reference_genome"
GRN_DIR="${PROJECT_DIR}/formatted_stability_GRNs_celltype"

# Where this run's formatted GRN ends up. Each method script hard-codes its own filename
# under ${GRN_DIR}/<Method>/, so the prefix has to be spelled out per method to know in
# advance whether a run is already done. DIRECTNET is deliberately absent: its script never
# writes a formatted GRN (there is no formatted_GRNs/DIRECTNET), so there is nothing to skip
# on and the check below is bypassed for it.
case "$METHOD" in
    CellOracle)  grn_prefix="celloracle" ;;
    LINGER)      grn_prefix="linger" ;;
    SCENIC_PLUS) grn_prefix="scenicplus" ;;
    FigR)        grn_prefix="figr" ;;
    Pando)       grn_prefix="pando" ;;
    *)           grn_prefix="" ;;
esac

if [ -n "${grn_prefix}" ]; then
    final_grn_file="${GRN_DIR}/${METHOD}/${grn_prefix}_${CELL_TYPE}_${RUN_SAMPLE_NAME}.tsv"
else
    final_grn_file=""
fi

# ===== FAILURE REPORTING =====
# Every task in the array appends one line to a shared failures file that sits next to the
# job's .out/.err files, so after the array finishes this one file lists exactly which runs
# died, which log to open, and which EXPERIMENT_LIST entries to re-run.
ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
slurm_log_dir="${PROJECT_DIR}/LOGS/run_tools_stability/run_tools_stability_${ARRAY_JOB_ID}"
failed_runs_file="${slurm_log_dir}/failed_runs_${ARRAY_JOB_ID}.txt"
mkdir -p "${slurm_log_dir}"

# ===== RESOURCE ALLOCATION =====
# Turns this job's #SBATCH header (-c / --mem, at the top of this file) into ALLOC_CPUS and
# ALLOC_MEM_*, which every method script consumes instead of choosing its own thread count.
# That is what makes the recorded resource usage comparable across methods.
source "${PROJECT_DIR}/src/common/resource_env.sh"
resource_env_summary

# ===== RESOURCE TRACKING =====
# The method scripts time their individual steps; this records what each method cost in
# total. Like the failures file, every task appends to one shared TSV next to the job's
# .out/.err files, so the whole array's resource usage lands in a single table.
source "${PROJECT_DIR}/src/common/resource_tracking.sh"
export RT_RECORD_FILE="${slurm_log_dir}/resource_usage_${ARRAY_JOB_ID}.tsv"
export RT_RAW_DIR="${slurm_log_dir}/resource_raw"
rt_init

failure_recorded=false
record_failure() {
    local reason="$1"
    # The EXIT trap fires after an explicit record_failure call, so guard against
    # logging the same task twice.
    if [ "${failure_recorded}" = "true" ]; then
        return 0
    fi
    failure_recorded=true

    # One printf of a short line is written atomically by the O_APPEND open, which keeps
    # the concurrent array tasks from interleaving their lines.
    printf 'task=%s\tmethod=%s\texperiment=%s\tsubsample=%s\treason=%s\tslurm_log=%s\tmethod_log=%s\n' \
        "${TASK_ID}" \
        "${METHOD}" \
        "${EXPERIMENT_CONFIG}" \
        "${SUBSAMPLE_NUM}" \
        "${reason}" \
        "${slurm_log_dir}/job_${TASK_ID}.err" \
        "${log_dir:-<not created>}" \
        >> "${failed_runs_file}"

    echo "Recorded failure in ${failed_runs_file}"
}

# ===== INTERMEDIATE RESULTS CLEANUP =====
# The formatted GRN is the only artifact worth keeping; this run's intermediates are large
# and disposable. That includes the filtered subsample RNA/ATAC inputs, which are written
# inside sample_result_dir for that reason. Failed runs get dropped too so they don't accumulate on disk -- the
# failures file above already records which log to open, and the method's own .log/.err
# under LOGS survive. Submit with KEEP_RESULTS=true
# (e.g. `sbatch --export=ALL,KEEP_RESULTS=true ...`) to keep them for debugging.
KEEP_RESULTS="${KEEP_RESULTS:-false}"
cleanup_results_dir() {
    # Guarded with :- because the trap is armed before sample_result_dir is assigned.
    if [ -z "${sample_result_dir:-}" ] || [ ! -d "${sample_result_dir}" ]; then
        return 0
    fi

    if [ "${KEEP_RESULTS}" = "true" ]; then
        echo "KEEP_RESULTS=true, leaving intermediate results dir: ${sample_result_dir}"
        return 0
    fi

    rm -rf "${sample_result_dir}"
    echo "Removed intermediate results dir: ${sample_result_dir}"

    # Prune the now-empty sample/cell_type parents; tasks still running keep their own dirs,
    # so this quietly no-ops for them. Bounded rmdir calls rather than `rmdir -p` so the
    # pruning can never walk up past the results root. (RESULTS_DIR has been reassigned to
    # the per-sample path by this point, so peel the parents off sample_result_dir instead.)
    # Three levels: the subsample, cell type and sample dirs.
    local prune_dir="${sample_result_dir}"
    for _ in 1 2 3; do
        prune_dir="${prune_dir%/*}"
        rmdir "${prune_dir}" 2>/dev/null || true
    done
}

on_exit() {
    local code=$?
    if [ "${code}" -ne 0 ]; then
        record_failure "exit_code=${code}"
    fi
    # In the trap so it runs however the script ends: set -e aborting on a failed method
    # never reaches the verification section below.
    cleanup_results_dir
}
trap on_exit EXIT

# ===== SKIP IF ALREADY DONE =====
# Reruns of the array shouldn't redo finished tasks. Submit with FORCE_RUN=true
# (e.g. `sbatch --export=ALL,FORCE_RUN=true ...`) to recompute regardless.
FORCE_RUN="${FORCE_RUN:-false}"
if [ "${FORCE_RUN}" != "true" ] && [ -n "${final_grn_file}" ] && [ -s "${final_grn_file}" ]; then
    echo "Formatted GRN already exists, skipping: ${final_grn_file}"
    echo "Set FORCE_RUN=true to rerun."
    exit 0
fi

mkdir -p "${RESULTS_DIR}" "${GRN_DIR}" "${DATA_DIR}"

# `|| true` because a no-match `ls` is a non-zero exit that set -e would otherwise turn
# into a bare exit 2, skipping the error message below.
rna_file=$(ls "${RAW_DATA_DIR}/${TISSUE}/${SAMPLE_NAME}/${CELL_TYPE}"/*RNA*.csv 2>/dev/null | head -n 1 || true)
atac_file=$(ls "${RAW_DATA_DIR}/${TISSUE}/${SAMPLE_NAME}/${CELL_TYPE}"/*ATAC*.csv 2>/dev/null | head -n 1 || true)

if [ -z "${rna_file}" ] || [ -z "${atac_file}" ]; then
    echo "ERROR: Could not find RNA or ATAC file"
    echo "RNA_FILE=${rna_file}"
    echo "ATAC_FILE=${atac_file}"
    record_failure "missing_input_files"
    exit 1
fi

# Keyed on METHOD as well as the sample: the array runs several methods for one sample
# concurrently, and cleanup_results_dir deletes this directory, so each task needs one of its
# own or a finishing job would wipe a running job's working files. (FigR, Pando and
# CellOracle re-append ${CELL_TYPE}/${SAMPLE_NAME} under whatever RESULTS_DIR they are given
# while LINGER, SCENIC_PLUS and DIRECTNET do not, so the depth below a method differs by
# method -- keying on METHOD keeps every one of them isolated regardless.)
# The subsample level keeps the 10 subsamples of one method/sample apart for the same reason.
sample_result_dir="${RESULTS_DIR}/${METHOD}/${SAMPLE_NAME}/${CELL_TYPE}/sub${SUBSAMPLE_NUM}"
mkdir -p "${sample_result_dir}"

log_dir="${PROJECT_DIR}/STABILITY_LOGS/${METHOD}/${SAMPLE_NAME}/${CELL_TYPE}/sub${SUBSAMPLE_NUM}"
mkdir -p "${log_dir}"

# ===== SUBSAMPLE INPUTS =====
# Each task writes its own copy of the full RNA/ATAC matrices, restricted to this subsample's
# barcodes. The copies go inside sample_result_dir, not /tmp: the node's /tmp is 20 GB and the
# larger ATAC subsamples are ~10 GB each. The EXIT trap removes them with the rest of the dir.
input_dir="${sample_result_dir}/subsample_input"
mkdir -p "${input_dir}"
sub_rna_file="${input_dir}/${CELL_TYPE}_RNA.csv"
sub_atac_file="${input_dir}/${CELL_TYPE}_ATAC.csv"

python "${PROJECT_DIR}/filter_rna_atac_stability.py" \
    --json_file "${STABILITY_JSON_FILE}" \
    --dataset "${TISSUE}" \
    --sample_name "${SAMPLE_NAME}" \
    --cell_type "${CELL_TYPE}" \
    --subsample_num "${SUBSAMPLE_NUM}" \
    --rna_file "${rna_file}" \
    --atac_file "${atac_file}" \
    --rna_out "${sub_rna_file}" \
    --atac_out "${sub_atac_file}" \
    || { record_failure "subsample_filter_failed"; exit 1; }

conda info --envs

# Shared environment
export PROJECT_DIR
export RAW_DATA_DIR
export DATA_DIR
export RESULTS_DIR="$sample_result_dir"
export REFERENCE_GENOME_DIR
export LOG_DIR="$log_dir"
export GRN_DIR
export CELL_TYPE
export SAMPLE_NAME="$RUN_SAMPLE_NAME"
export SUBSAMPLE_NUM
export SPECIES
export RNA_FILE="$sub_rna_file"
export ATAC_FILE="$sub_atac_file"

# The case only resolves which script to run; rt_run_method below does the running, so the
# resource tracking wraps every method identically instead of being repeated per branch.
case "$METHOD" in
    CellOracle)  method_script="${PROJECT_DIR}/src/Celloracle/run_CellOracle.sh" ;;
    DIRECTNET)   method_script="${PROJECT_DIR}/src/DIRECTNET/run_DIRECTNET.sh" ;;
    LINGER)      method_script="${PROJECT_DIR}/src/LINGER/run_linger.sh" ;;
    SCENIC_PLUS) method_script="${PROJECT_DIR}/src/SCENIC_PLUS/run_scenic_plus.sh" ;;
    FigR)        method_script="${PROJECT_DIR}/src/FigR/run_FigR.sh" ;;
    Pando)       method_script="${PROJECT_DIR}/src/Pando/run_pando.sh" ;;
    *)
        echo "ERROR: Unknown method: ${METHOD}"
        exit 1
        ;;
esac

# Runs the method under /usr/bin/time -v and appends its totals to RT_RECORD_FILE, returning
# the method's own exit status so set -e still aborts on a failed method.
rt_run_method "$METHOD" "$method_script" \
    "${log_dir}/${METHOD}.log" \
    "${log_dir}/${METHOD}.err"

echo ""
echo "===== RESOURCE USAGE (${METHOD} total) ====="
rt_summary_last
echo "  record    : ${RT_RECORD_FILE}"
echo "  raw       : ${RT_LAST_RAW}"

# ===== VERIFY FORMATTED GRN =====
# The methods write their formatted GRN themselves, so there is nothing to move here -- but a
# method can still exit 0 without producing one, and that run needs re-running just as much as
# a crashed one. Skipped for DIRECTNET, which produces no formatted GRN by design.
echo ""
if [ -n "${final_grn_file}" ]; then
    if [ -s "${final_grn_file}" ]; then
        echo "Formatted GRN: ${final_grn_file}"
        # sample_result_dir is removed by the EXIT trap, on this path and on failures alike.
    else
        echo "WARNING: expected a formatted GRN at ${final_grn_file}, found none"
        record_failure "no_formatted_grn"
    fi
else
    echo "No formatted GRN expected for ${METHOD}, skipping verification"
fi

echo "Completed ${METHOD} for ${SAMPLE_NAME}/${CELL_TYPE} (task ${TASK_ID})"
