# multiGRNtools

## Tool Installation
### LINGER Installation

Install Conda environment
```bash
conda env create -f linger_environment.yml
```

---

### SCENIC+ Installation

Handled during setup automatically. To manually set up the Conda environment, run:
```bash
conda env create -f scenicplus_environment.yml
```
---

### DIRECTNET Installation

Install Conda environment
```bash
conda env create -f directnet_environment.yml
source activate directnet_env
```

Install R dependencies
```bash
module load R/4.3.2
Rscript install_directnet.R
```

---

### CellOracle Installation

Install Conda environment
```bash
conda env create -f celloracle_environment.yml
```

