#!/bin/bash -l
#SBATCH --job-name="install_directnet"
#SBATCH --output=LOGS/install_directnet_%j.log
#SBATCH --error=LOGS/install_directnet_%j.log
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=16G

source activate directnet_env

Rscript src/DIRECTNET/install_directnet.R