library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install xfun version 0.55 to avoid the error "Error in getNamespaceExports("xfun") : object 'attr' not found"
remotes::install_version(
  "xfun",
  version = "0.55",
  lib = "/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3",
  repos = "https://cloud.r-project.org",
  upgrade = "never"
)

# Install htmlTable (required by Hmisc)
install.packages("htmlTable", lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install Hmisc (required by biovizBase)
install.packages("Hmisc", lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install biovizBase (required by Gviz)
BiocManager::install(c("biovizBase"), lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install Gviz, GenomicRanges, and rtracklayer (required by cicero)
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"), lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install cicero (required by DIRECT-NET)
devtools::install_github("cole-trapnell-lab/cicero-release", ref="monocle3", lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install Signac (required by DIRECT-NET)
devtools::install_github("timoast/signac", ref = "develop", lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

# Install DIRECT-NET
devtools::install_github("zhanglhbioinfor/DIRECT-NET", lib="/gpfs/Home/esm5360/R/x86_64-conda-linux-gnu-library/4.3")

library(cicero)
library(Signac)
library(DIRECTNET)