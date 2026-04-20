#!/bin/bash -l

#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=23
#SBATCH --mem-per-cpu=16G

set -euo pipefail

# ==========================================
#             USER VARIABLES
# ==========================================

# Conda environment name
CONDA_ENV_NAME="LINGER"

METHOD="scNN"

ACTIVEF='ReLU'

BASE_DIR="${PROJECT_DIR}/src/LINGER"
LINGER_RESULTS_DIR="${RESULTS_DIR}/LINGER"
TSS_MOTIF_INFO_PATH="${DATA_DIR}/LINGER_OTHER_SPECIES_TF_MOTIF_DATA/provide_data/"

if [[ "$SPECIES" == "human" ]]; then
    GENOME="hg38"
elif [[ "$SPECIES" == "mouse" ]]; then
    GENOME="mm10"
fi

# ==========================================
#             SETUP FUNCTIONS
# ==========================================

validate_critical_variables() {
    echo "[INFO] Validating required variables..."
    local required_vars=(
        SAMPLE_NAME
        RNA_FILE
        ATAC_FILE
        BASE_DIR
        DATA_DIR
        RESULTS_DIR
    )
    for var in "${required_vars[@]}"; do
        if [ -z "${!var:-}" ]; then
            echo "[ERROR] Required variable '$var' is not set."
            exit 1
        fi
    done
}

check_for_running_jobs() {
    echo ""
    echo "[INFO] Checking for running jobs with the same name..."
    if [ -n "${SLURM_JOB_NAME:-}" ]; then
        local running_jobs
        running_jobs=$(squeue --name="$SLURM_JOB_NAME" --noheader | wc -l)
        if [ "$running_jobs" -gt 1 ]; then
            echo "[ERROR] Another job with the same name is already running. Exiting."
            exit 1
        else
            echo "    - No conflicting jobs detected."
        fi
    else
        echo "    - Not running under SLURM, skipping job check."
    fi
}

check_tools() {
    echo ""
    echo "[INFO] Checking for required tools (python3, conda)..."
    for tool in python3 conda; do
        if ! command -v "$tool" &> /dev/null; then
            echo "[ERROR] Required tool '$tool' is not installed or not in PATH."
            exit 1
        else
            echo "    - Found $tool"
        fi
    done
}

activate_conda_env() {
    echo ""
    echo "[INFO] Activating Conda environment '$CONDA_ENV_NAME'..."
    if ! conda activate "$CONDA_ENV_NAME"; then
        echo "[ERROR] Could not activate Conda environment '$CONDA_ENV_NAME'."
        exit 1
    fi
    echo "    - Successfully activated Conda environment."
}

check_input_files() {
    echo ""
    echo "[INFO] Checking if input RNA and ATAC files exist..."
    for file in "$RNA_FILE" "$ATAC_FILE"; do
        if [ ! -f "$file" ]; then
            echo "[ERROR] Missing input file: $file"
            exit 1
        else
            echo "    - Found $file"
        fi
    done
}

setup_directories() {
    echo ""
    echo "[INFO] Setting up necessary directories..."
    mkdir -p "$LINGER_RESULTS_DIR"
    echo "    - Directories created."
}

set_slurm_job_name() {
    echo ""
    echo "[INFO] Setting dynamic SLURM job name..."
    scontrol update JobID="$SLURM_JOB_ID" JobName="LINGER_${SAMPLE_NAME}"
}

check_or_install_tss_locations() {
    # Checks for the TSS location information for species other than hg38
    echo "  1) Checking for other species TSS location directory..."

    if [ ! -d "${TSS_MOTIF_INFO_PATH}" ]; then
        echo "[WARN] TSS location not found at ${TSS_MOTIF_INFO_PATH}"
        echo "    - Attempting to download and extract..."

        mkdir -p "${DATA_DIR}/LINGER_OTHER_SPECIES_TF_MOTIF_DATA"

        confirm_code=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate \
            'https://drive.usercontent.google.com/download?id=1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV' -O- | \
            sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p')

        wget --load-cookies /tmp/cookies.txt \
            "https://drive.usercontent.google.com/download?export=download&confirm=${confirm_code}&id=1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV" \
            -O "${DATA_DIR}/LINGER_OTHER_SPECIES_TF_MOTIF_DATA/provide_data.tar.gz"

        rm -rf /tmp/cookies.txt

        tar -xzf "${DATA_DIR}/LINGER_OTHER_SPECIES_TF_MOTIF_DATA/provide_data.tar.gz" -C "${DATA_DIR}/LINGER_OTHER_SPECIES_TF_MOTIF_DATA/"
        echo "    - Downloaded and extracted TSS location information to ${TSS_MOTIF_INFO_PATH}"
    else
        echo "    - TSS location information exists"
    fi

    # Download hg38 TSS locations from ENSEMBL biomart using the provided script if not present
    hg38_tss_file="${TSS_MOTIF_INFO_PATH}/TSS_hg38.txt"
    if [ ! -f "${hg38_tss_file}" ]; then
        python3 "${BASE_DIR}/download_hg38_ensembl_tss.py" \
            --save_file "${hg38_tss_file}"
        
        if [ -f "${hg38_tss_file}" ]; then
            echo "    - Successfully downloaded hg38 TSS location file."


        else
            echo "[ERROR] Failed to download hg38 TSS location file."
            exit 1
        fi
    fi

    # Add hg38 TSS information to the genome_map_homer.txt file used by the Homer motif finding step
    HOMER_GENOME_MAP_FILE="${TSS_MOTIF_INFO_PATH}/genome_map_homer.txt"
    touch "${HOMER_GENOME_MAP_FILE}"
    if ! grep -Eq '^hg38([[:space:]]|$)' "${HOMER_GENOME_MAP_FILE}"; then
        printf 'hg38\thuman\tHomo_sapiens\tGRCh38\t102\tMammal\n' >> "${HOMER_GENOME_MAP_FILE}"
        echo "    - Added hg38 information to genome_map_homer.txt"
    fi
}


check_homer() {
    echo "  2) Checking Homer installation..."

    # Check if 'homer' is installed
    if ! command -v homer &> /dev/null; then
        echo "[ERROR] Homer is not installed or not in PATH."
        echo "    - Installing Homer..."
        
        # Homer installation directory
        HOMER_INSTALL_DIR="${DATA_DIR}/homer"

        mkdir -p "$HOMER_INSTALL_DIR"
        cd "$HOMER_INSTALL_DIR" || exit 1

        # Download and install Homer
        wget http://homer.ucsd.edu/homer/configureHomer.pl
        perl configureHomer.pl -install

        # Add to PATH temporarily (you should add permanently too)
        export PATH="$HOMER_INSTALL_DIR/bin:$PATH"

        echo "    - Homer installed at $HOMER_INSTALL_DIR"
    else
        echo "    - Homer already installed."
    fi

    HOMER_BIN=$(dirname "$(which homer)")
    HOMER_BASE=$(realpath "${HOMER_BIN}/..")
    HOMER_CONFIG="${HOMER_BASE}/share/homer/config.txt"

    # Check if Homer is installed
    if [ ! -f "$HOMER_CONFIG" ]; then
        echo "[ERROR] HOMER config.txt not found at $HOMER_CONFIG"
        echo "Make sure Homer is properly installed."
        exit 1
    fi

    # Check if GENOME exists in config.txt
    if grep -q "${GENOME}" "$HOMER_CONFIG"; then
        echo "[INFO] Homer genome '${GENOME}' already installed."
    else
        echo "[WARN] Homer genome '${GENOME}' not found. Installing..."
        perl "${HOMER_BASE}/share/homer/configureHomer.pl" -install "${GENOME}" \
            2> "${LOG_DIR}/install_homer_species.log"
        
        if grep -q "${GENOME}" "$HOMER_CONFIG"; then
            echo "[INFO] Successfully installed '${GENOME}'."
        else
            echo "[ERROR] Failed to install '${GENOME}'. Check log: ${LOG_DIR}/install_homer_species.log"
            exit 1
        fi
    fi

    # Add homer/bin to PATH if not already
    if [[ ":$PATH:" != *":${HOMER_BASE}/homer/bin:"* ]]; then
        echo "[INFO] Adding Homer bin to PATH."
        export PATH="${PATH}:${HOMER_BASE}/homer/bin"
    fi
}

determine_num_cpus() {
    echo ""
    echo "[INFO] Checking the number of CPUs available for parallel processing"
    if [ -z "${SLURM_CPUS_PER_TASK:-}" ]; then
        if command -v nproc &> /dev/null; then
            TOTAL_CPUS=$(nproc --all)
            case $TOTAL_CPUS in
                [1-15]) IGNORED_CPUS=1 ;;  # Reserve 1 CPU for <=15 cores
                [16-31]) IGNORED_CPUS=2 ;; # Reserve 2 CPUs for <=31 cores
                *) IGNORED_CPUS=4 ;;       # Reserve 4 CPUs for >=32 cores
            esac
            NUM_CPU=$((TOTAL_CPUS - IGNORED_CPUS))
            echo "    - Running locally. Detected $TOTAL_CPUS CPUs, reserving $IGNORED_CPUS for system tasks. Using $NUM_CPU CPUs."
        else
            NUM_CPU=1  # Fallback
            echo "    - Running locally. Unable to detect CPUs, defaulting to $NUM_CPU CPU."
        fi
    else
        NUM_CPU=${SLURM_CPUS_PER_TASK}
        echo "    - Running on SLURM. Number of CPUs allocated: ${NUM_CPU}"
    fi
}


# ==========================================
#             PIPELINE STEPS
# ==========================================

# Function to run each Python step with timing and memory tracking
run_step() {
  step_name=$1
  script_path=$2
  shift 2  # Remove the first two arguments
  echo ""
  echo "Running $step_name..."
  /usr/bin/time -v python3 "$script_path" "$@" 2> "${LOG_DIR}/${step_name}_time_mem.log"
  echo "    Done!"
}

run_pipeline() {
    run_step "Step_010.Linger_Load_Data" "${BASE_DIR}/Step_010.Linger_Load_Data.py" \
        --rna_data_path "$RNA_FILE" \
        --atac_data_path "$ATAC_FILE" \
        --data_dir "$DATA_DIR" \
        --sample_data_dir "$LINGER_RESULTS_DIR" \
        --organism "$SPECIES" \
        --cell_type "$CELL_TYPE" \
        --genome "$GENOME" \
        --method "$METHOD"

    run_step "Step_020.Linger_Training" "${BASE_DIR}/Step_020.Linger_Training.py" \
        --tss_motif_info_path "$TSS_MOTIF_INFO_PATH" \
        --genome "$GENOME" \
        --method "$METHOD" \
        --sample_data_dir "$LINGER_RESULTS_DIR" \
        --activef "$ACTIVEF" \
        --organism "$SPECIES" 

    run_step "Step_030.Create_Cell_Population_GRN" "${BASE_DIR}/Step_030.Create_Cell_Population_GRN.py" \
        --tss_motif_info_path "$TSS_MOTIF_INFO_PATH" \
        --genome "$GENOME" \
        --method "$METHOD" \
        --sample_data_dir "$LINGER_RESULTS_DIR" \
        --activef "$ACTIVEF" \
        --organism "$SPECIES" \
        --num_cpu "$NUM_CPU"

    run_step "Step_040.Homer_Motif_Finding" "${BASE_DIR}/Step_040.Homer_Motif_Finding.py" \
        --tss_motif_info_path "$TSS_MOTIF_INFO_PATH" \
        --sample_data_dir "$LINGER_RESULTS_DIR" \
        --genome "$GENOME"

    run_step "Step_050.Create_Cell_Type_GRN" "${BASE_DIR}/Step_050.Create_Cell_Type_GRN.py" \
        --tss_motif_info_path "$TSS_MOTIF_INFO_PATH" \
        --genome "$GENOME" \
        --method "$METHOD" \
        --sample_data_dir "$LINGER_RESULTS_DIR" \
        --cell_type "$CELL_TYPE" \
        --organism "$SPECIES" \
        --num_cpu "$NUM_CPU"

    GRN_FILE="${LINGER_RESULTS_DIR}/cell_type_specific_trans_regulatory_${CELL_TYPE}.txt"
    run_step "Step_060.Format_Inferred_GRN" "${BASE_DIR}/Step_060.Format_Inferred_Grn.py" \
        --output_dir "${GRN_DIR}/LINGER/" \
        --inferred_grn_file "${GRN_FILE}" \
        --cell_type "$CELL_TYPE" \
        --sample_name "$SAMPLE_NAME"
}

# ==========================================
#               MAIN
# ==========================================

echo "===== RUNNING VALIDATION CHECKS ====="
validate_critical_variables
check_for_running_jobs
check_tools
check_input_files
activate_conda_env
setup_directories
set_slurm_job_name
determine_num_cpus
check_or_install_tss_locations
check_homer

echo ""
echo "===== CHECKS COMPLETE: STARTING MAIN PIPELINE ====="
run_pipeline
