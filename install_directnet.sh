#!/bin/bash -l
#SBATCH --job-name="install_directnet"
#SBATCH --output=LOGS/DIRECTNET/install_directnet.log
#SBATCH --error=LOGS/DIRECTNET/install_directnet.log
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=16G

conda activate directnet_env

CONDA_ENV_LIB="$CONDA_PREFIX/lib/R/library"

Rscript src/DIRECTNET/install_directnet.R \
    "$CONDA_ENV_LIB"

