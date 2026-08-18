#!/bin/bash -l
#SBATCH --job-name=CellOracle
#SBATCH --output=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/CellOracle/CellOracle_%A.txt
#SBATCH --error=/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/multiGRNtools/LOGS/CellOracle/CellOracle_%A.err
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G

GENOME="${GENOME:-Human}"    # Default; may be overridden by SPECIES below

## Map SPECIES from submit script to CellOracle genome mode
if [[ -n "$SPECIES" ]]; then
  case "${SPECIES,,}" in
    human) GENOME="Human" ;;
    mouse) GENOME="Mouse" ;;
    *)
      echo "ERROR: SPECIES must be 'human' or 'mouse' when provided. Got: '$SPECIES'"
      exit 1
      ;;
  esac
fi

## ── Conda + modules ───────────────────────────────────────────
source activate celloracle_env
module load bedtools

if [[ -z "$CONDA_PREFIX" ]]; then
  echo "ERROR: conda environment activation failed (CONDA_PREFIX is empty)"
  exit 1
fi

PYTHON_BIN="$CONDA_PREFIX/bin/python"
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "ERROR: Python executable not found in conda env: $PYTHON_BIN"
  exit 1
fi

## Provide newer conda libstdc++ for user R packages compiled with newer ABI.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
if [[ -f "$CONDA_PREFIX/lib/libstdc++.so.6" ]]; then
  export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
else
  echo "ERROR: Expected libstdc++ not found at $CONDA_PREFIX/lib/libstdc++.so.6"
  exit 1
fi

## ── Base paths ────────────────────────────────────────────────
BASE_DIR="${PROJECT_DIR}/src/Celloracle"
SCRIPT_DIR=$BASE_DIR
CELL_ORACLE_RESULTS_DIR="${RESULTS_DIR}/${CELL_TYPE}/${SAMPLE_NAME}/CellOracle"
# Honour the LOG_DIR the parent runner exports, so stability runs land in
# STABILITY_LOGS/.../subsample_N instead of all subsamples sharing one file.
# Falls back to the standalone default when invoked directly via sbatch.
LOG_DIR="${LOG_DIR:-${PROJECT_DIR}/LOGS/CellOracle/${CELL_TYPE}/${SAMPLE_NAME}}"

mkdir -p "$LOG_DIR"

## ── Genome-specific settings ──────────────────────────────────
if [[ "$GENOME" == "Human" ]]; then
  GENOME_REF="hg38"
  CHROM_SIZES="${SCRIPT_DIR}/hg38.fa.sizes"

elif [[ "$GENOME" == "Mouse" ]]; then
  GENOME_REF="mm10"
  CHROM_SIZES="${SCRIPT_DIR}/mm10.chrom.sizes"

else
  echo "ERROR: GENOME must be 'Human' or 'Mouse'. Got: '$GENOME'"
  exit 1
fi

## ── Validate required paths ───────────────────────────────────
if [[ -z "$CELL_TYPE" || -z "$SAMPLE_NAME" ]]; then
  echo "ERROR: CELL_TYPE and SAMPLE_NAME must be provided via --export"
  exit 1
fi
if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "ERROR: Chromosome sizes file not found: $CHROM_SIZES"
  exit 1
fi
if [[ ! -d "$REFERENCE_GENOME_DIR/$GENOME_REF" ]]; then
  echo "ERROR: Genome directory not found: $REFERENCE_GENOME_DIR/$GENOME_REF"
  exit 1
fi

## genomepy/gimmemotifs background generation may still use genomepy default path.
## Mirror shared genome directories into that default location via symlinks.
GENOMEPY_DEFAULT_DIR="$HOME/.local/share/genomes"
mkdir -p "$GENOMEPY_DEFAULT_DIR"

if [[ -d "$REFERENCE_GENOME_DIR/hg38" ]]; then
  ln -sfn "$REFERENCE_GENOME_DIR/hg38" "$GENOMEPY_DEFAULT_DIR/hg38"
fi
if [[ -d "$REFERENCE_GENOME_DIR/mm10" ]]; then
  ln -sfn "$REFERENCE_GENOME_DIR/mm10" "$GENOMEPY_DEFAULT_DIR/mm10"
fi

if [[ ! -f "$GENOMEPY_DEFAULT_DIR/$GENOME_REF/$GENOME_REF.fa" && ! -f "$GENOMEPY_DEFAULT_DIR/$GENOME_REF/$GENOME_REF.fa.gz" ]]; then
  echo "ERROR: Genome FASTA not found for $GENOME_REF under $GENOMEPY_DEFAULT_DIR/$GENOME_REF"
  exit 1
fi

if [[ -z "$RNA_FILE" || -z "$ATAC_FILE" ]]; then
  echo "ERROR: RNA_FILE and ATAC_FILE are required via --export"
  exit 1
fi
if [[ ! -f "$RNA_FILE" ]]; then
  echo "ERROR: RNA_FILE not found: $RNA_FILE"
  exit 1
fi
if [[ ! -f "$ATAC_FILE" ]]; then
  echo "ERROR: ATAC_FILE not found: $ATAC_FILE"
  exit 1
fi

SAMPLE_DIR="$(dirname "$RNA_FILE")"

echo "========================================"
echo "  Genome       : $GENOME"
echo "  Genome ref   : $GENOME_REF"
echo "  Chrom sizes  : $CHROM_SIZES"
echo "  Script dir   : $SCRIPT_DIR"
echo "  Rscript      : $(command -v Rscript)"
echo "  Python       : $PYTHON_BIN"
echo "  LD_PRELOAD   : ${LD_PRELOAD:-<empty>}"
echo "  Cell type    : $CELL_TYPE"
echo "  Sample name  : $SAMPLE_NAME"
echo "  Sample dir   : $SAMPLE_DIR"
echo "  RNA_FILE     : $RNA_FILE"
echo "  ATAC_FILE    : $ATAC_FILE"
echo "  Genomes dir  : $REFERENCE_GENOME_DIR"
echo "  Genomepy dir : $GENOMEPY_DEFAULT_DIR"
echo "  Results dir  : $CELL_ORACLE_RESULTS_DIR"
echo "========================================"
SAMPLE="$SAMPLE_NAME"

# ALLOC_CPUS is set by src/common/resource_env.sh from the runner's #SBATCH header, so every
# method gets the identical budget. The fallback keeps this script runnable standalone.
# Re-exported because CellOracle's parallelism lives inside step7_final_integration.py rather
# than in a command-line argument -- that step read a hard-coded n_jobs=4 before this, so
# CellOracle was benchmarked on 4 cores while every other method used the full allocation.
NUM_THREADS="${ALLOC_CPUS:-${SLURM_CPUS_PER_TASK:-1}}"
export ALLOC_CPUS="${NUM_THREADS}"

echo "========================================"
echo "  Sample       : $SAMPLE"
echo "  Sample dir   : $SAMPLE_DIR"
echo "  Node         : $SLURMD_NODENAME"
echo "  Num threads  : $NUM_THREADS"
echo "  Start        : $(date)"
echo "========================================"

## ══════════════════════════════════════════════════════════════
##  run_condition — execute the full 4-step CellOracle pipeline
##
##  Arguments:
##    $1 rna_file   : full path to RNA CSV
##    $2 atac_file  : full path to ATAC CSV
##    $3 out_dir    : full output directory for this sample
##    $4 cell_label : label written as cell_type in AnnData
##                    (optional; defaults to $SAMPLE)
## ══════════════════════════════════════════════════════════════
run_condition() {
  local rna_file="$1"
  local atac_file="$2"
  local out_dir="$3"
  local grn_dir="$4"
  local cell_label="${5:-$SAMPLE}"

  echo ""
  echo "     RNA        : $rna_file"
  echo "     ATAC       : $atac_file"
  echo "     Cell label : $cell_label"
  echo "     Out dir    : $out_dir"

  if [[ ! -f "$rna_file"  ]]; then echo "  SKIP: RNA file not found : $rna_file";  return; fi
  if [[ ! -f "$atac_file" ]]; then echo "  SKIP: ATAC file not found: $atac_file"; return; fi

  mkdir -p "$out_dir"

  ## Intermediate file paths — all anchored to out_dir
  local out_atac_prefix="${out_dir}/atac"
  local out_rna_prefix="${out_dir}/rna"
  local all_peak="${out_atac_prefix}_all_peaks.csv"
  local cicero_connection="${out_atac_prefix}_cicero_connections.csv"
  local int_atac_input="${out_atac_prefix}_base_GRN.parquet"
  local int_rna_input="${out_rna_prefix}_processed.h5ad"

  {
    echo "=== CellOracle ==="
    echo "  RNA        : $rna_file"
    echo "  ATAC       : $atac_file"
    echo "  Cell label : $cell_label"
    echo "  Out dir    : $out_dir"
    echo "  Genome ref : $GENOME_REF"
    echo "  Start      : $(date)"

    ## [1/4] Preprocess ATAC → all_peaks.csv + cicero_connections.csv
    echo ""
    echo "  [1/4] step1_preprocess.R ..."
    /usr/bin/time -v Rscript "${SCRIPT_DIR}/step1_preprocess.R" \
      "$atac_file" \
      "$out_atac_prefix" \
      "$CHROM_SIZES" \
      || { echo "  FAILED at step1_preprocess.R"; exit 1; }

    ## [2/4] TSS annotation + motif scanning → base GRN parquet
    echo ""
    echo "  [2/4] step3_TSS_annot.py ..."
    /usr/bin/time -v "$PYTHON_BIN" "${SCRIPT_DIR}/step3_TSS_annot.py" \
      "$all_peak" \
      "$cicero_connection" \
      "$GENOME_REF" \
      "$out_atac_prefix" \
      "$REFERENCE_GENOME_DIR" \
      || { echo "  FAILED at step3_TSS_annot.py"; exit 1; }

    ## [3/4] Process RNA → processed .h5ad
    echo ""
    echo "  [3/4] steps_all_RNA.py ..."
    /usr/bin/time -v "$PYTHON_BIN" "${SCRIPT_DIR}/steps_all_RNA.py" \
      "$rna_file" \
      "$out_rna_prefix" \
      "$cell_label" \
      || { echo "  FAILED at steps_all_RNA.py"; exit 1; }

    ## [4/4] Integrate RNA + ATAC → GRN per cluster
    echo ""
    echo "  [4/4] step7_final_integration.py ..."
    /usr/bin/time -v "$PYTHON_BIN" "${SCRIPT_DIR}/step7_final_integration.py" \
      "$int_rna_input" \
      "$int_atac_input" \
      "$out_dir" \
      || { echo "  FAILED at step7_final_integration.py"; exit 1; }

    echo ""
    echo "  Done: $(date)"

    echo ""
    echo "Formatting GRN..."
    /usr/bin/time -v "$PYTHON_BIN" "${SCRIPT_DIR}/format_celloracle_grn.py" \
      "$out_dir/output_GRN/celloracle_network_${SAMPLE}_final_GRN.csv" \
      "$grn_dir/CellOracle/celloracle_${CELL_TYPE}_${SAMPLE}.tsv"

  } > "$LOG_DIR/CellOracle_run.log" 2>&1 \
    && echo "  Done" \
    || echo "  FAILED : see $LOG_DIR/CellOracle_run.log"
}

run_condition \
  "$RNA_FILE" \
  "$ATAC_FILE" \
  "$CELL_ORACLE_RESULTS_DIR" \
  "$GRN_DIR" \
  "$SAMPLE"
  

echo ""
echo "========================================"
echo "  Sample $SAMPLE complete : $(date)"
echo "========================================"
