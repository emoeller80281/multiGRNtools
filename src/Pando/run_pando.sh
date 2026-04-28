#!/bin/bash -l
#SBATCH --job-name=Pando
#SBATCH --output=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/Pando/Pando_%A.txt
#SBATCH --error=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/Pando/Pando_%A.err
#SBATCH --time=12:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

GENOME="${GENOME:-Human}"

if [[ -n "${SPECIES:-}" ]]; then
  case "${SPECIES,,}" in
    human) GENOME="Human" ;;
    mouse) GENOME="Mouse" ;;
    *)
      echo "ERROR: SPECIES must be 'human' or 'mouse' when provided. Got: '$SPECIES'"
      exit 1
      ;;
  esac
fi

## ── modules ───────────────────────────────────────────
source activate pando_env

if [[ -z "$CONDA_PREFIX" ]]; then
  echo "ERROR: conda environment activation failed (CONDA_PREFIX is empty)"
  exit 1
fi

PYTHON_BIN="$CONDA_PREFIX/bin/python"
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python executable not found in conda env: $PYTHON_BIN"
  exit 1
fi

# Provide newer conda libstdc++ for user R packages compiled with newer ABI.
if [[ -z "$CONDA_PREFIX" ]]; then
  echo "ERROR: conda environment activation failed (CONDA_PREFIX is empty)"
  exit 1
fi
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

if [[ ! -f "$CONDA_PREFIX/lib/libstdc++.so.6" ]]; then
  echo "ERROR: Expected libstdc++ not found at $CONDA_PREFIX/lib/libstdc++.so.6"
  exit 1
fi

## ── Base paths ────────────────────────────────────────────────
BASE_DIR="${PROJECT_DIR}/src/Pando"
SCRIPT_DIR=$BASE_DIR
Pando_RESULTS_DIR="${RESULTS_DIR}/${CELL_TYPE}/${SAMPLE_NAME}/Pando"
LOG_DIR="${PROJECT_DIR}/LOGS/Pando/${CELL_TYPE}/${SAMPLE_NAME}"

mkdir -p "$LOG_DIR" "$Pando_RESULTS_DIR"


## Provide newer conda libstdc++ for user R packages compiled with newer ABI.

if [[ -z "${CELL_TYPE:-}" || -z "${SAMPLE_NAME:-}" ]]; then
  echo "ERROR: CELL_TYPE and SAMPLE_NAME must be provided via --export"
  exit 1
fi
if [[ -z "${RNA_FILE:-}" || -z "${ATAC_FILE:-}" ]]; then
  echo "ERROR: RNA_FILE and ATAC_FILE are required via --export"
  exit 1
fi
# if [[ -z "${REFERENCE_GENOME_DIR:-}" ]]; then
#   echo "ERROR: REFERENCE_GENOME_DIR must be provided via --export"
#   exit 1
# fi
if [[ ! -f "$RNA_FILE" ]]; then
  echo "ERROR: RNA_FILE not found: $RNA_FILE"
  exit 1
fi
if [[ ! -f "$ATAC_FILE" ]]; then
  echo "ERROR: ATAC_FILE not found: $ATAC_FILE"
  exit 1
fi


## ── Genome-specific settings ──────────────────────────────────
if [[ "$GENOME" == "Human" ]]; then
  GENOME_REF="hg38"
elif [[ "$GENOME" == "Mouse" ]]; then
  GENOME_REF="mm10"
else
  echo "ERROR: GENOME must be 'Human' or 'Mouse'. Got: '$GENOME'"
  exit 1
fi

SAMPLE_DIR="$(dirname "$RNA_FILE")"
SAMPLE="$SAMPLE_NAME"
NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "========================================"
echo "  Genome      : $GENOME"
echo "  Genome ref  : $GENOME_REF"
echo "  Project dir  : $PROJECT_DIR"
echo "  Results dir  : $RESULTS_DIR"
echo "  Script dir   : $SCRIPT_DIR"
echo "  Log dir      : $LOG_DIR"
echo "  Cell type    : $CELL_TYPE"
echo "  Sample name  : $SAMPLE_NAME"
echo "  Sample dir   : $SAMPLE_DIR"
echo "  RNA_FILE     : $RNA_FILE"
echo "  ATAC_FILE    : $ATAC_FILE"
echo "  Num threads   : $NUM_THREADS"
echo "========================================"


run_pando() {
  if [[ ! -f "$RNA_FILE" ]]; then echo "  SKIP: RNA file not found : $RNA_FILE"; return; fi
  if [[ ! -f "$ATAC_FILE" ]]; then echo "  SKIP: ATAC file not found: $ATAC_FILE"; return; fi

  mkdir -p "$Pando_RESULTS_DIR"

  {
    echo "=== Starting Pando ==="
    echo "  Start     : $(date)"

    /usr/bin/time -v Rscript "${SCRIPT_DIR}/Pando.R" \
      "$RNA_FILE" \
      "$ATAC_FILE" \
      "$Pando_RESULTS_DIR" \
      "$SAMPLE" \
      "$GENOME_REF" \
      "$SCRIPT_DIR" \
      "$NUM_THREADS" \
      || { echo "  FAILED at Pando.R"; exit 1; }

    echo ""
    echo "  Done: $(date)"
  } > "$LOG_DIR/Pando_run.log" 2>&1 \
    && echo "  Done"
}

run_pando

echo ""
echo "Formatting inferred GRN"
PANDO_RAW_GRN_CSV="${Pando_RESULTS_DIR}/${SAMPLE_NAME}_filtered_network.csv"
PANDO_FORMATTED_GRN_TSV="${GRN_DIR}/Pando/pando_${CELL_TYPE}_${SAMPLE_NAME}.tsv"

mkdir -p "$(dirname "$PANDO_FORMATTED_GRN_TSV")"

python "${SCRIPT_DIR}/format_pando_grn.py" \
  "$PANDO_RAW_GRN_CSV" \
  "$PANDO_FORMATTED_GRN_TSV"

echo ""
echo "========================================"
echo "  Sample $SAMPLE complete : $(date)"
echo "========================================"
