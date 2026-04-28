#!/bin/bash -l
#SBATCH --job-name=FigR
#SBATCH --output=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/FigR/FigR_%A.txt
#SBATCH --error=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/FigR/FigR_%A.err
#SBATCH --time=12:00:00
#SBATCH -p memory
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
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
source activate figr_env

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
BASE_DIR="${PROJECT_DIR}/src/FigR"
SCRIPT_DIR=$BASE_DIR
FigR_RESULTS_DIR="${RESULTS_DIR}/${CELL_TYPE}/${SAMPLE_NAME}/FigR"
LOG_DIR="${PROJECT_DIR}/LOGS/FigR/${CELL_TYPE}/${SAMPLE_NAME}"

mkdir -p "$LOG_DIR" "$FigR_RESULTS_DIR"


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

# GENOME_DIR="${REFERENCE_GENOME_DIR}/${GENOME_REF}"
# TSS_FILE="${GENOME_DIR}/gene_tss.bed"

# if [[ "$GENOME_REF" == "hg38" ]]; then
#   GTF_FILE="${GENOME_DIR}/Homo_sapiens.GRCh38.113.gtf"
# else
#   GTF_FILE="${GENOME_DIR}/Mus_musculus.GRCm39.113.gtf"
# fi

# if [[ ! -d "$GENOME_DIR" ]]; then
#   echo "ERROR: Genome directory not found: $GENOME_DIR"
#   exit 1
# fi
# if [[ ! -f "$TSS_FILE" ]]; then
#   echo "ERROR: TSS file not found: $TSS_FILE"
#   exit 1
# fi
# if [[ ! -f "$GTF_FILE" ]]; then
#   echo "ERROR: GTF file not found: $GTF_FILE"
#   exit 1
# fi

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
# echo "  Genome dir   : $GENOME_DIR"
# echo "  TSS file     : $TSS_FILE"
# echo "  GTF file     : $GTF_FILE"
echo "  Num threads   : $NUM_THREADS"
echo "========================================"


run_figr() {
  if [[ ! -f "$RNA_FILE" ]]; then echo "  SKIP: RNA file not found : $RNA_FILE"; return; fi
  if [[ ! -f "$ATAC_FILE" ]]; then echo "  SKIP: ATAC file not found: $ATAC_FILE"; return; fi

  mkdir -p "$FigR_RESULTS_DIR"

  {
    echo "=== Starting FigR ==="
    echo "  Start     : $(date)"

    /usr/bin/time -v Rscript "${SCRIPT_DIR}/FigR.R" \
      "$RNA_FILE" \
      "$ATAC_FILE" \
      "$FigR_RESULTS_DIR" \
      "$SAMPLE" \
      "$GENOME_REF" \
      "$NUM_THREADS" \
      || { echo "  FAILED at FigR.R"; exit 1; }

    echo ""
    echo "  Done: $(date)"
  } > "$LOG_DIR/FigR_run.log" 2>&1 \
    && echo "  Done"
}

run_figr

echo ""
echo "========================================"
echo "  Sample $SAMPLE complete : $(date)"
echo "========================================"
