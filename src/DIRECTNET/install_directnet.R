personal_lib <- "/gpfs/Home/esm5360/miniconda3/envs/directnet_env/lib/R/library"

remotes::install_github(
  "cole-trapnell-lab/monocle3",
  lib = personal_lib,
  upgrade = "never"
)

remotes::install_github(
  "cole-trapnell-lab/cicero-release",
  ref = "monocle3",
  lib = personal_lib,
  upgrade = "never"
)

remotes::install_github(
  "zhanglhbioinfor/DIRECT-NET",
  lib = "/gpfs/Home/esm5360/miniconda3/envs/directnet_env/lib/R/library",
  upgrade = "never"
)

# Rscript -e 'remotes::install_github("zhanglhbioinfor/DIRECT-NET", upgrade="never", lib = "/gpfs/Home/esm5360/miniconda3/envs/directnet_env/lib/R/library")'

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# # Install xfun version 0.55 to avoid the error "Error in getNamespaceExports("xfun") : object 'attr' not found"
# remotes::install_version(
#   "xfun",
#   version = "0.55",
#   lib = personal_lib,
#   repos = "https://cloud.r-project.org",
#   upgrade = "never"
# )

# # Install htmlTable (required by Hmisc)
# install.packages("htmlTable", lib=personal_lib)

# # Install Hmisc (required by biovizBase)
# install.packages("Hmisc", lib=personal_lib)

# # Install biovizBase (required by Gviz)
# BiocManager::install(c("biovizBase"), lib=personal_lib)

# # Install Gviz, GenomicRanges, and rtracklayer (required by cicero)
# BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"), lib=personal_lib)

# # Install grr (required by monocle3)
# remotes::install_github("cran/grr", lib=personal_lib)

# # Install monocle3 (required by cicero)
# devtools::install_github('cole-trapnell-lab/monocle3', lib=personal_lib)

# # Install cicero (required by DIRECT-NET)
# devtools::install_github("cole-trapnell-lab/cicero-release", ref="monocle3", lib=personal_lib)

# # Install Signac (required by DIRECT-NET)
# devtools::install_github("timoast/signac", ref = "develop", lib=personal_lib)

# install.packages("xgboost", lib=personal_lib)

# # Install DIRECT-NET
# devtools::install_github("zhanglhbioinfor/DIRECT-NET", lib=personal_lib)

# # ---

# library(cicero)
# library(Signac)
# library(DIRECTNET)