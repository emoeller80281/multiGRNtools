#!/bin/bash -l
#SBATCH --job-name="multiGRNtools_celltype"
#SBATCH --output="LOGS/run_tools/run_tools_%A/job_%a.out"
#SBATCH --error="LOGS/run_tools/run_tools_%A/job_%a.err"
#SBATCH --time=24:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 12
#SBATCH --mem=128G
#SBATCH --exclude=psh01com1ochm[01-03]
#SBATCH --array=0-219%20

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
    # "ENCODE_colon|colon_adult_ENCDO856ZOJ|Colonic_epithelium|human"     #  5,666 cells
    # "ENCODE_colon|colon_adult_ENCDO856ZOJ|Smooth_muscle_cells|human"    #    529 cells
    # "ENCODE_colon|colon_adult_ENCDO856ZOJ|Fibroblasts|human"            #    266 cells
    # "ENCODE_colon|colon_adult_ENCDO856ZOJ|Macrophages|human"            #    123 cells
    # "ENCODE_colon|colon_adult_ENCDO856ZOJ|Endothelial_cells|human"      #    121 cells
    # "ENCODE_adrenal|adrenal_adult_ENCDO221IPH|Adrenal_cortex_cells|human"   #  2,640 cells
    # "ENCODE_adrenal|adrenal_adult_ENCDO221IPH|Fibroblasts|human"            #    506 cells
    # "ENCODE_adrenal|adrenal_adult_ENCDO221IPH|Endothelial_cells|human"      #    492 cells
    # "ENCODE_adrenal|adrenal_adult_ENCDO221IPH|Macrophages|human"            #    192 cells
    # "ENCODE_adrenal|adrenal_adult_ENCDO221IPH|Smooth_muscle_cells|human"    #     51 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Cholangiocytes|human"        #  1,956 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Fibroblasts|human"           #  1,180 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Macrophages|human"           #    477 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Schwann_cells|human"         #    359 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Endothelial_cells|human"     #    334 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|T_cells|human"               #    200 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Smooth_muscle_cells|human"   #    131 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Mast_cells|human"            #     80 cells
    # "ENCODE_bile_duct|bile_duct_adult_ENCDO575WHY|Hepatocytes|human"           #     73 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|Pancreatic_exocrine_cells|human"    #  9,610 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|Fibroblasts|human"                  #    309 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|Pancreatic_endocrine_cells|human"   #    219 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|Endothelial_cells|human"            #    181 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|T_cells|human"                      #     22 cells
    # "ENCODE_pancreas|pancreas_adult_ENCDO221IPH|Macrophages|human"                  #     18 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO411EVD|Fallopian_tube_epithelium|human"   #  4,084 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO411EVD|Fibroblasts|human"                 #    345 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO411EVD|Smooth_muscle_cells|human"         #    150 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO411EVD|Endothelial_cells|human"           #    132 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO411EVD|Macrophages|human"                 #     68 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|Fibroblasts|human"         #  2,597 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|Endothelial_cells|human"   #  1,374 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|Cardiomyocytes|human"      #  1,296 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|Smooth_muscle_cells|human" #    998 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|Macrophages|human"         #    845 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO575WHY|T_cells|human"             #    180 cells
    # "ENCODE_lung|lung_adult_ENCDO528BHB|Fibroblasts|human"            #    635 cells
    # "ENCODE_lung|lung_adult_ENCDO528BHB|Endothelial_cells|human"      #    496 cells
    # "ENCODE_lung|lung_adult_ENCDO528BHB|Macrophages|human"            #    127 cells
    # "ENCODE_lung|lung_adult_ENCDO528BHB|Airway_epithelium|human"      #     78 cells
    # "ENCODE_lung|lung_adult_ENCDO528BHB|Smooth_muscle_cells|human"    #     58 cells
    # "GSE277326_retina_macula_rpe|Retina_donor1|Rod_photoreceptors|human"    #  2,017 cells
    # "GSE277326_retina_macula_rpe|Retina_donor1|Cone_photoreceptors|human"   #     54 cells
    # "GSE277326_retina_macula_rpe|Retina_donor1|Astrocytes|human"            #     29 cells
    # "GSE277326_retina_macula_rpe|Retina_donor1|Microglia|human"             #     10 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Skeletal_myocytes|human"     #  1,968 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Fibroblasts|human"           #    259 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Satellite_cells|human"       #     90 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Macrophages|human"           #     47 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Endothelial_cells|human"     #     40 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr1|Smooth_muscle_cells|human"   #     21 cells
    # "GSE207335_epidermis|S1|Keratinocytes|human"      #  1,733 cells
    # "GSE207335_epidermis|S1|Langerhans_cells|human"   #     31 cells
    # "GSE207335_epidermis|S1|T_cells|human"            #     23 cells
    # "GSE207335_epidermis|S1|Melanocytes|human"        #     15 cells
    # "GSE194122_bone_marrow|s1d1|T_cells|human"           #  1,728 cells
    # "GSE194122_bone_marrow|s1d1|Erythroid_cells|human"   #    726 cells
    # "GSE194122_bone_marrow|s1d1|Monocytes|human"         #    725 cells
    # "GSE194122_bone_marrow|s1d1|NK_cells|human"          #    523 cells
    # "GSE194122_bone_marrow|s1d1|B_cells|human"           #    500 cells
    # "GSE194122_bone_marrow|s1d1|HSPCs|human"             #    170 cells
    # "GSE194122_bone_marrow|s1d1|Dendritic_cells|human"   #    156 cells
    # "GSE194122_bone_marrow|s1d1|pDCs|human"              #     59 cells
    # "GSE194122_bone_marrow|s1d1|Plasma_cells|human"      #     51 cells
    # "GSE319044_lung_spleen_immune|Spleen-SMO-04|B_cells|human"   #  1,229 cells
    # "GSE319044_lung_spleen_immune|Spleen-SMO-04|T_cells|human"   #    790 cells
    "ENCODE_adrenal|adrenal_adult_ENCDO856ZOJ|Adrenal_cortex_cells|human"  #  6,341 cells
    "ENCODE_adrenal|adrenal_adult_ENCDO856ZOJ|Macrophages|human"           #    584 cells
    "ENCODE_adrenal|adrenal_adult_ENCDO856ZOJ|Fibroblasts|human"           #    348 cells
    "ENCODE_adrenal|adrenal_adult_ENCDO856ZOJ|Endothelial_cells|human"     #    340 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Colonic_epithelium|human"   #  2,178 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Smooth_muscle_cells|human"  #    233 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Plasma_cells|human"         #    224 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Fibroblasts|human"          #    206 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|T_cells|human"              #    161 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Macrophages|human"          #    144 cells
    "ENCODE_colon|colon_adult_ENCDO186XRB|Endothelial_cells|human"    #    141 cells
    "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|Fallopian_tube_epithelium|human"  #  3,761 cells
    "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|Fibroblasts|human"                #    267 cells
    "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|Endothelial_cells|human"          #    174 cells
    "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|Smooth_muscle_cells|human"        #    163 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|Macrophages|human"                #     38 cells
    # "ENCODE_fallopian_tube|fallopian_tube_adult_ENCDO856ZOJ|T_cells|human"                    #      9 cells
    "ENCODE_heart|heart_lv_adult_ENCDO411EVD|Cardiomyocytes|human"       #  5,359 cells
    "ENCODE_heart|heart_lv_adult_ENCDO411EVD|Fibroblasts|human"          #  2,414 cells
    "ENCODE_heart|heart_lv_adult_ENCDO411EVD|Macrophages|human"          #  1,512 cells
    "ENCODE_heart|heart_lv_adult_ENCDO411EVD|Smooth_muscle_cells|human"  #  1,228 cells
    "ENCODE_heart|heart_lv_adult_ENCDO411EVD|Endothelial_cells|human"    #  1,138 cells
    # "ENCODE_heart|heart_lv_adult_ENCDO411EVD|T_cells|human"              #     84 cells
    "ENCODE_lung|lung_child_ENCDO575EGL|Fibroblasts|human"          #    891 cells
    "ENCODE_lung|lung_child_ENCDO575EGL|Endothelial_cells|human"    #    316 cells
    "ENCODE_lung|lung_child_ENCDO575EGL|Macrophages|human"          #    222 cells
    # "ENCODE_lung|lung_child_ENCDO575EGL|Smooth_muscle_cells|human"  #     99 cells
    # "ENCODE_lung|lung_child_ENCDO575EGL|Plasma_cells|human"         #     36 cells
    "ENCODE_pancreas|pancreas_child_ENCDO575EGL|Pancreatic_exocrine_cells|human"  #  4,260 cells
    "ENCODE_pancreas|pancreas_child_ENCDO575EGL|Macrophages|human"                #    583 cells
    "ENCODE_pancreas|pancreas_child_ENCDO575EGL|Fibroblasts|human"                #    387 cells
    "ENCODE_pancreas|pancreas_child_ENCDO575EGL|T_cells|human"                    #    100 cells
    # "ENCODE_pancreas|pancreas_child_ENCDO575EGL|Endothelial_cells|human"          #     66 cells
    # "ENCODE_pancreas|pancreas_child_ENCDO575EGL|Smooth_muscle_cells|human"        #     42 cells
    "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|Fibroblasts|human"          #  1,336 cells
    "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|Smooth_muscle_cells|human"  #  1,295 cells
    "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|Macrophages|human"          #    488 cells
    "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|Endothelial_cells|human"    #    358 cells
    "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|T_cells|human"              #    107 cells
    # "ENCODE_vena_cava|vena_cava_adult_ENCDO221IPH|Schwann_cells|human"        #     25 cells
    "GSE194122_bone_marrow|s1d2|T_cells|human"          #  3,146 cells
    "GSE194122_bone_marrow|s1d2|NK_cells|human"         #  1,089 cells
    "GSE194122_bone_marrow|s1d2|Monocytes|human"        #    768 cells
    "GSE194122_bone_marrow|s1d2|pDCs|human"             #    241 cells
    "GSE194122_bone_marrow|s1d2|B_cells|human"          #    229 cells
    "GSE194122_bone_marrow|s1d2|Erythroid_cells|human"  #    149 cells
    "GSE207335_epidermis|S2|Keratinocytes|human"     #  2,177 cells
    # "GSE207335_epidermis|S2|Langerhans_cells|human"  #     69 cells
    "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Skeletal_myocytes|human"    #  1,892 cells
    "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Fibroblasts|human"          #    241 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Endothelial_cells|human"    #     70 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Macrophages|human"          #     64 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Satellite_cells|human"      #     50 cells
    # "GSE240061_muscle_pre_exercise|Ex-1_Pre_tr2|Smooth_muscle_cells|human"  #     21 cells
    "GSE277326_retina_macula_rpe|Retina_donor2|Rod_photoreceptors|human"   #  2,801 cells
    # "GSE277326_retina_macula_rpe|Retina_donor2|Cone_photoreceptors|human"  #     79 cells
    "GSE319044_lung_spleen_immune|Spleen-SMO-1|B_cells|human"  #  4,819 cells
    "GSE319044_lung_spleen_immune|Spleen-SMO-1|T_cells|human"  #    163 cells
)

NUM_EXPERIMENTS=${#EXPERIMENT_LIST[@]}
NUM_METHODS=${#METHOD_LIST[@]}
TOTAL_TASKS=$((NUM_EXPERIMENTS * NUM_METHODS))

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ "$TASK_ID" -ge "$TOTAL_TASKS" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${TASK_ID} exceeds TOTAL_TASKS=${TOTAL_TASKS}"
    exit 1
fi

EXPERIMENT_ID=$((TASK_ID / NUM_METHODS))
METHOD_ID=$((TASK_ID % NUM_METHODS))

EXPERIMENT_CONFIG="${EXPERIMENT_LIST[$EXPERIMENT_ID]}"
METHOD="${METHOD_LIST[$METHOD_ID]}"

IFS='|' read -r TISSUE SAMPLE_NAME CELL_TYPE SPECIES <<< "$EXPERIMENT_CONFIG"

echo "TASK_ID=${TASK_ID}"
echo "EXPERIMENT_ID=${EXPERIMENT_ID}"
echo "METHOD_ID=${METHOD_ID}"
echo "METHOD=${METHOD}"
echo "TISSUE=${TISSUE}"
echo "SAMPLE_NAME=${SAMPLE_NAME}"
echo "CELL_TYPE=${CELL_TYPE}"
echo "SPECIES=${SPECIES}"


# ===== PATH CONFIGURATION =====
PROJECT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
RAW_DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/MUON_FILTERED_COUNT_DATASETS/hg38"
DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/data"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
REFERENCE_GENOME_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.SINGLE_CELL_GRN_INFERENCE.MOELLER/data/genome_data/reference_genome"
GRN_DIR="${PROJECT_DIR}/hg38_formatted_GRNs_celltype"

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
    final_grn_file="${GRN_DIR}/${METHOD}/${grn_prefix}_${CELL_TYPE}_${SAMPLE_NAME}.tsv"
else
    final_grn_file=""
fi

# ===== FAILURE REPORTING =====
# Every task in the array appends one line to a shared failures file that sits next to the
# job's .out/.err files, so after the array finishes this one file lists exactly which runs
# died, which log to open, and which EXPERIMENT_LIST entries to re-run.
ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
slurm_log_dir="${PROJECT_DIR}/LOGS/run_tools/run_tools_${ARRAY_JOB_ID}"
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
    printf 'task=%s\tmethod=%s\texperiment=%s\treason=%s\tslurm_log=%s\tmethod_log=%s\n' \
        "${TASK_ID}" \
        "${METHOD}" \
        "${EXPERIMENT_CONFIG}" \
        "${reason}" \
        "${slurm_log_dir}/job_${TASK_ID}.err" \
        "${log_dir:-<not created>}" \
        >> "${failed_runs_file}"

    echo "Recorded failure in ${failed_runs_file}"
}

# ===== INTERMEDIATE RESULTS CLEANUP =====
# The formatted GRN is the only artifact worth keeping; this run's intermediates are large
# and disposable. Failed runs get dropped too so they don't accumulate on disk -- the
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
    # Two levels here rather than the stability script's three: there is no subsample dir.
    local prune_dir="${sample_result_dir}"
    for _ in 1 2; do
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
sample_result_dir="${RESULTS_DIR}/${METHOD}/${SAMPLE_NAME}/${CELL_TYPE}"
mkdir -p "${sample_result_dir}"

log_dir="${PROJECT_DIR}/LOGS/${METHOD}/${SAMPLE_NAME}/${CELL_TYPE}"
mkdir -p "${log_dir}"

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
export SAMPLE_NAME
export SPECIES
export RNA_FILE="$rna_file"
export ATAC_FILE="$atac_file"

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
