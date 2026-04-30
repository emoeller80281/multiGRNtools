#!/bin/bash -l
#SBATCH --job-name="multiGRNtools_2d"
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=48:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH --array=0-100%5

set -euo pipefail

# ===== SELECTED METHODS =====
METHOD_LIST=(
    "Pando"
    "FigR"
    # "SCENIC_PLUS"
    # "LINGER"
    # "DIRECTNET"
    # "CellOracle"
)

# ===== SAMPLE CONFIGURATION =====
EXPERIMENT_LIST=(
    # "mESC|E7.5_rep1|mouse|mESC"
    # "mESC|E7.5_rep2|mouse|mESC"
    # "mESC|E8.5_rep1|mouse|mESC"
    # "mESC|E8.5_rep2|mouse|mESC"
    # "Macrophage|buffer_1|human|Macrophage"
    "Macrophage|buffer_2|human|Macrophage"
    # "Macrophage|buffer_3|human|Macrophage"
    # "Macrophage|buffer_4|human|Macrophage"
    # "K562|sample_1|human|K562"
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

IFS='|' read -r CELL_TYPE SAMPLE_NAME SPECIES RAW_CELL_TYPE <<< "$EXPERIMENT_CONFIG"

echo "TASK_ID=${TASK_ID}"
echo "EXPERIMENT_ID=${EXPERIMENT_ID}"
echo "METHOD_ID=${METHOD_ID}"
echo "METHOD=${METHOD}"
echo "CELL_TYPE=${CELL_TYPE}"
echo "SAMPLE_NAME=${SAMPLE_NAME}"
echo "SPECIES=${SPECIES}"
echo "RAW_CELL_TYPE=${RAW_CELL_TYPE}"

# ===== PATH CONFIGURATION =====
PROJECT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
RAW_DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/MUON_FILTERED_COUNT_DATASETS"
DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/data"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
REFERENCE_GENOME_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.SINGLE_CELL_GRN_INFERENCE.MOELLER/data/genome_data/reference_genome"
GRN_DIR="${PROJECT_DIR}/formatted_GRNs"

mkdir -p "${RESULTS_DIR}" "${GRN_DIR}" "${DATA_DIR}"

rna_file=$(ls "${RAW_DATA_DIR}/${RAW_CELL_TYPE}/${SAMPLE_NAME}"/*RNA*.csv 2>/dev/null | head -n 1)
atac_file=$(ls "${RAW_DATA_DIR}/${RAW_CELL_TYPE}/${SAMPLE_NAME}"/*ATAC*.csv 2>/dev/null | head -n 1)

if [ -z "${rna_file}" ] || [ -z "${atac_file}" ]; then
    echo "ERROR: Could not find RNA or ATAC file"
    echo "RNA_FILE=${rna_file}"
    echo "ATAC_FILE=${atac_file}"
    exit 1
fi

sample_result_dir="${RESULTS_DIR}/${CELL_TYPE}/${SAMPLE_NAME}"
mkdir -p "${sample_result_dir}"

log_dir="${PROJECT_DIR}/LOGS/${METHOD}/${CELL_TYPE}/${SAMPLE_NAME}"
mkdir -p "${log_dir}"

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

case "$METHOD" in
    CellOracle)
        bash "${PROJECT_DIR}/src/Celloracle/run_CellOracle.sh" \
            > "${log_dir}/CellOracle.log" \
            2> "${log_dir}/CellOracle.err"
        ;;

    DIRECTNET)
        bash "${PROJECT_DIR}/src/DIRECTNET/run_DIRECTNET.sh" \
            > "${log_dir}/DIRECTNET.log" \
            2> "${log_dir}/DIRECTNET.err"
        ;;

    LINGER)
        bash "${PROJECT_DIR}/src/LINGER/run_linger.sh" \
            > "${log_dir}/LINGER.log" \
            2> "${log_dir}/LINGER.err"
        ;;

    SCENIC_PLUS)
        bash "${PROJECT_DIR}/src/SCENIC_PLUS/run_scenic_plus.sh" \
            > "${log_dir}/SCENIC_PLUS.log" \
            2> "${log_dir}/SCENIC_PLUS.err"
        ;;

    FigR)
        bash "${PROJECT_DIR}/src/FigR/run_FigR.sh" \
            > "${log_dir}/FigR.log" \
            2> "${log_dir}/FigR.err"
        ;;

    Pando)
        bash "${PROJECT_DIR}/src/Pando/run_pando.sh" \
            > "${log_dir}/Pando.log" \
            2> "${log_dir}/Pando.err"
        ;;

    *)
        echo "ERROR: Unknown method: ${METHOD}"
        exit 1
        ;;
esac

echo "Completed ${METHOD} for ${CELL_TYPE}/${SAMPLE_NAME}"