#!/bin/bash -l
#SBATCH --job-name="install_directnet"
#SBATCH --output=LOGS/install_directnet_%j.log
#SBATCH --error=LOGS/install_directnet_%j.log
#SBATCH --time=08:00:00
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=16G

# conda activate directnet_env
module load rstudio

## Provide newer conda libstdc++ for user R packages compiled with newer ABI.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
if [[ -f "$CONDA_PREFIX/lib/libstdc++.so.6" ]]; then
  export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
else
  echo "ERROR: Expected libstdc++ not found at $CONDA_PREFIX/lib/libstdc++.so.6"
  exit 1
fi

# Remove stale lock files
rm -rf ~/R/x86_64-conda-linux-gnu-library/4.3/00LOCK*

Rscript src/DIRECTNET/install_directnet.R