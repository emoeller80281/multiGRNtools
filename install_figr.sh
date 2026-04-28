#!/bin/bash -l
#SBATCH --job-name="install_figr"
#SBATCH --output=LOGS/FigR/install_figr.log
#SBATCH --error=LOGS/FigR/install_figr.log
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=16G

conda activate figr_env

CONDA_ENV_LIB="$CONDA_PREFIX/lib/R/library"

Rscript src/FigR/install_figr.R \
    "$CONDA_ENV_LIB"


