options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  download.file.method = "libcurl"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# ---- CRAN ----
install.packages(c(
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
))

# ---- Bioconductor ----
BiocManager::install(c(
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "chromVAR",
  "motifmatchr"
))

# ---- Genome packages (choose based on your runs) ----
BiocManager::install(c(
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10"
))

# ---- Cicero (required by DIRECTNET) ----
remotes::install_github(
  "cole-trapnell-lab/cicero-release",
  ref = "monocle3"
)

# ---- DIRECTNET ----
remotes::install_github("zhanglhbioinfor/DIRECT-NET")