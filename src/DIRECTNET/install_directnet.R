# ----------------------------
# User library setup
# ----------------------------
dir.create("~/R/library", recursive = TRUE, showWarnings = FALSE)
.libPaths(c("~/R/library", .libPaths()))

options(download.file.method = "libcurl")

# ----------------------------
# Install/load BiocManager first
# ----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Now that BiocManager exists, point repos to Bioconductor + CRAN
options(repos = BiocManager::repositories())

cat("R version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
print(BiocManager::repositories())

# ----------------------------
# Pin to the R-4.3-compatible Bioconductor release
# ----------------------------
BiocManager::install(version = "3.18", ask = FALSE)
options(repos = BiocManager::repositories())

# ----------------------------
# CRAN helpers first
# ----------------------------
install.packages(
  c("remotes", "fs"),
  repos = "https://cloud.r-project.org"
)

# ----------------------------
# Core Bioconductor deps needed by DIRECTNET
# Install these before DIRECTNET
# ----------------------------
BiocManager::install(
  c(
    "sparseMatrixStats",
    "GenomicRanges",
    "GenomeInfoDb",
    "rtracklayer",
    "chromVAR",
    "motifmatchr"
  ),
  ask = FALSE,
  update = FALSE
)

# ----------------------------
# CRAN packages
# ----------------------------
install.packages(
  c(
    "Seurat",
    "Signac",
    "patchwork",
    "dplyr",
    "ggplot2",
    "Matrix",
    "FNN",
    "slam",
    "xgboost",
    "stringr",
    "RColorBrewer"
  ),
  repos = "https://cloud.r-project.org"
)

# ----------------------------
# Genome annotation packages
# Keep both only if you truly need both species
# ----------------------------
BiocManager::install(
  c(
    "EnsDb.Hsapiens.v86",
    "BSgenome.Hsapiens.UCSC.hg38",
    "EnsDb.Mmusculus.v79",
    "BSgenome.Mmusculus.UCSC.mm10"
  ),
  ask = FALSE,
  update = FALSE
)

# ----------------------------
# Cicero (required by DIRECTNET)
# ----------------------------
remotes::install_github(
  "cole-trapnell-lab/cicero-release",
  ref = "monocle3",
  upgrade = "never",
  dependencies = FALSE
)

# ----------------------------
# DIRECTNET
# Only after chromVAR + motifmatchr are installed
# ----------------------------
remotes::install_github(
  "zhanglhbioinfor/DIRECT-NET",
  upgrade = "never",
  dependencies = FALSE
)