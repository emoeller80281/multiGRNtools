#!/bin/bash -l
#SBATCH --job-name="submit_multiple_scmultipredict_jobs"
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --array=0%10

PROJECT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
RAW_DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/MUON_FILTERED_COUNT_DATASETS"
DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/data"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
REFERENCE_GENOME_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.SINGLE_CELL_GRN_INFERENCE.MOELLER/data/genome_data/reference_genome"
GRN_DIR="${PROJECT_DIR}/formatted_GRNs"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${GRN_DIR}"
mkdir -p "${DATA_DIR}"

RUN_CELLORACLE=false
RUN_LINGER=true

EXPERIMENT_LIST=(
    # "mESC|E7.5_rep1|mouse|mESC"
    # "mESC|E7.5_rep2|mouse|mESC"
    # "mESC|E8.5_rep1|mouse|mESC"
    # "mESC|E8.5_rep2|mouse|mESC"

    # "Macrophage|buffer_1|human|Macrophage"
    # "Macrophage|buffer_2|human|Macrophage"
    # "Macrophage|buffer_3|human|Macrophage"
    # "Macrophage|buffer_4|human|Macrophage"

    # "iPSC|WT_D13_rep1|human|iPSC"

    "K562|sample_1|human|K562"
)

# ==========================================
#        TASK SELECTION
# ==========================================
# Get the current experiment based on SLURM_ARRAY_TASK_ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ ${TASK_ID} -ge ${#EXPERIMENT_LIST[@]} ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID (${TASK_ID}) exceeds number of experiments (${#EXPERIMENT_LIST[@]})"
    exit 1
fi

EXPERIMENT_CONFIG="${EXPERIMENT_LIST[$TASK_ID]}"

IFS='|' read -r CELL_TYPE SAMPLE_NAME SPECIES RAW_CELL_TYPE <<< "$EXPERIMENT_CONFIG"

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

rna_file="${RAW_DATA_DIR}/${RAW_CELL_TYPE}/${SAMPLE_NAME}/${SAMPLE_NAME}_RNA.csv"
atac_file="${RAW_DATA_DIR}/${RAW_CELL_TYPE}/${SAMPLE_NAME}/${SAMPLE_NAME}_ATAC.csv"

sample_result_dir="${RESULTS_DIR}/${CELL_TYPE}/${SAMPLE_NAME}"
mkdir -p "${sample_result_dir}"

if [ "$RUN_CELLORACLE" = true ]; then
    echo "Submitting CellOracle job for ${CELL_TYPE} - ${SAMPLE_NAME} (Task ID: ${ARRAY_TASK_ID})"

    log_dir="${PROJECT_DIR}/LOGS/CellOracle/${CELL_TYPE}/${SAMPLE_NAME}"
    mkdir -p "${log_dir}"

    sbatch \
        --export=PROJECT_DIR="$PROJECT_DIR",RAW_DATA_DIR="$RAW_DATA_DIR",RESULTS_DIR="$sample_result_dir",REFERENCE_GENOME_DIR="$REFERENCE_GENOME_DIR",CELL_TYPE="$CELL_TYPE",SAMPLE_NAME="$SAMPLE_NAME",SPECIES="$SPECIES",RNA_FILE="$rna_file",ATAC_FILE="$atac_file" \
        --job-name="SCMULTI_PREDICT_CELLORACLE_${CELL_TYPE}_${SAMPLE_NAME}" \
        --output=${log_dir}/CellOracle.log \
        --error=${log_dir}/CellOracle.err \
        "${PROJECT_DIR}/src/Celloracle/run_CellOracle.sh"

if [ "$RUN_LINGER" = true ]; then
    echo "Submitting LINGER job for ${CELL_TYPE} - ${SAMPLE_NAME} (Task ID: ${ARRAY_TASK_ID})"

    log_dir="${PROJECT_DIR}/LOGS/LINGER/${CELL_TYPE}/${SAMPLE_NAME}"
    mkdir -p "${log_dir}"

    sbatch \
        --export=PROJECT_DIR="$PROJECT_DIR",DATA_DIR="$DATA_DIR",RAW_DATA_DIR="$RAW_DATA_DIR",RESULTS_DIR="$sample_result_dir",LOG_DIR="$log_dir",GRN_DIR="$GRN_DIR",CELL_TYPE="$CELL_TYPE",SAMPLE_NAME="$SAMPLE_NAME",SPECIES="$SPECIES",RNA_FILE="$rna_file",ATAC_FILE="$atac_file" \
        --job-name="SCMULTI_PREDICT_LINGER_${CELL_TYPE}_${SAMPLE_NAME}" \
        --output=${log_dir}/LINGER.log \
        --error=${log_dir}/LINGER.err \
        "${PROJECT_DIR}/src/LINGER/run_linger.sh"