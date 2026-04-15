#!/bin/bash -l
#SBATCH --job-name="submit_multiple_scmultipredict_jobs"
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --array=0%10

PROJECT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"
RAW_DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/MUON_FILTERED_COUNT_DATASETS"
RESULTS_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools"


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

# Run CellOracle for each sample
sbatch \
    --export=PROJECT_DIR="$PROJECT_DIR",RAW_DATA_DIR="$RAW_DATA_DIR",RESULTS_DIR="$RESULTS_DIR",CELL_TYPE="$CELL_TYPE",SAMPLE_NAME="$SAMPLE_NAME",SPECIES="$SPECIES",RNA_FILE="$rna_file",ATAC_FILE="$atac_file" \
    --job-name="SCMULTI_PREDICT_${CELL_TYPE}_${SAMPLE_NAME}" \
    --output=${PROJECT_DIR}/LOGS/CellOracle/${CELL_TYPE}/${SAMPLE_NAME}/CellOracle_${ARRAY_JOB_ID}_${ARRAY_TASK_ID}.log \
    --error=${PROJECT_DIR}/LOGS/CellOracle/${CELL_TYPE}/${SAMPLE_NAME}/CellOracle_${ARRAY_JOB_ID}_${ARRAY_TASK_ID}.err \
    "${PROJECT_DIR}/Celloracle/Mouse/run_CellOracle1.slurm"
