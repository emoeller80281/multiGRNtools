#!/bin/bash -l
#SBATCH --job-name="install_pando"
#SBATCH --output=LOGS/Pando/install_pando.log
#SBATCH --error=LOGS/Pando/install_pando.log
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=16G

conda activate pando_env

CONDA_ENV_LIB="$CONDA_PREFIX/lib/R/library"

Rscript src/Pando/install_pando.R \
    "$CONDA_ENV_LIB"


