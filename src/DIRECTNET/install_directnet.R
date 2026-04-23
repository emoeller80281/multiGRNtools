# ----------------------------
# User library setup
# ----------------------------
user_lib <- Sys.getenv("R_LIBS_USER")

if (!nzchar(user_lib)) {
  user_lib <- "~/R/%p-library/%v"
  user_lib <- path.expand(
    gsub(
      "%p", R.version$platform,
      gsub("%v", paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = "."),
           user_lib)
    )
  )
} else {
  user_lib <- path.expand(user_lib)
}

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, setdiff(.libPaths(), user_lib)))

cat("Using R library:", user_lib, "\n")
print(.libPaths())

# ----------------------------
# Force C++14 where needed
# ----------------------------
dir.create("~/.R", recursive = TRUE, showWarnings = FALSE)
makevars_path <- path.expand("~/.R/Makevars")

makevars_lines <- c(
  "CXX11STD = -std=gnu++14",
  "CXX14STD = -std=gnu++14",
  "CXX17STD = -std=gnu++17"
)

if (file.exists(makevars_path)) {
  existing <- readLines(makevars_path, warn = FALSE)
  keep <- existing[!grepl("^(CXX11STD|CXX14STD|CXX17STD)\\s*=", existing)]
  writeLines(c(keep, makevars_lines), makevars_path)
} else {
  writeLines(makevars_lines, makevars_path)
}

cat("Using Makevars at:", makevars_path, "\n")

# ----------------------------
# Helpers
# ----------------------------
is_installed <- function(pkg, lib = .libPaths()[1]) {
  requireNamespace(pkg, quietly = TRUE, lib.loc = lib)
}

check_all_installed <- function(pkgs, lib = .libPaths()[1]) {
  vapply(pkgs, is_installed, logical(1), lib = lib)
}

install_if_missing_cran <- function(pkgs, repos = "https://cloud.r-project.org", lib = .libPaths()[1], ...) {
  missing <- pkgs[!check_all_installed(pkgs, lib = lib)]
  if (length(missing) == 0) {
    cat("All CRAN packages already installed:", paste(pkgs, collapse = ", "), "\n")
    return(invisible(TRUE))
  }

  cat("Installing missing CRAN packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = repos, lib = lib, ...)

  still_missing <- missing[!check_all_installed(missing, lib = lib)]
  if (length(still_missing) > 0) {
    stop("These CRAN packages failed to install: ", paste(still_missing, collapse = ", "))
  }

  invisible(TRUE)
}

install_if_missing_bioc <- function(pkgs, lib = .libPaths()[1], ..., ask = FALSE, update = FALSE, dependencies = TRUE) {
  missing <- pkgs[!check_all_installed(pkgs, lib = lib)]
  if (length(missing) == 0) {
    cat("All Bioconductor packages already installed:", paste(pkgs, collapse = ", "), "\n")
    return(invisible(TRUE))
  }

  cat("Installing missing Bioconductor packages:", paste(missing, collapse = ", "), "\n")
  BiocManager::install(
    missing,
    ask = ask,
    update = update,
    lib = lib,
    dependencies = dependencies,
    ...
  )

  still_missing <- missing[!check_all_installed(missing, lib = lib)]
  if (length(still_missing) > 0) {
    stop("These Bioconductor packages failed to install: ", paste(still_missing, collapse = ", "))
  }

  invisible(TRUE)
}

install_github_if_missing <- function(pkg_name, repo, ref = NULL, force = FALSE, lib = .libPaths()[1], ...) {
  if (!force && is_installed(pkg_name, lib = lib)) {
    cat("GitHub package already installed:", pkg_name, "\n")
    return(invisible(TRUE))
  }

  cat("Installing GitHub package:", pkg_name, "from", repo, "\n")

  if (is.null(ref)) {
    remotes::install_github(repo, force = force, lib = lib, ...)
  } else {
    remotes::install_github(repo, ref = ref, force = force, lib = lib, ...)
  }

  if (!is_installed(pkg_name, lib = lib)) {
    stop("GitHub package failed to install: ", pkg_name)
  }

  invisible(TRUE)
}

# Useful when a package had a broken prior install
reinstall_bioc <- function(pkgs, lib = .libPaths()[1]) {
  cat("Reinstalling Bioconductor packages:", paste(pkgs, collapse = ", "), "\n")
  BiocManager::install(
    pkgs,
    ask = FALSE,
    update = FALSE,
    force = TRUE,
    lib = lib,
    dependencies = TRUE
  )
}

reinstall_cran <- function(pkgs, repos = "https://cloud.r-project.org", lib = .libPaths()[1]) {
  cat("Reinstalling CRAN packages:", paste(pkgs, collapse = ", "), "\n")
  install.packages(pkgs, repos = repos, lib = lib)
}

# ----------------------------
# Install BiocManager first if missing
# ----------------------------
if (!is_installed("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = .libPaths()[1])
}

options(download.file.method = "libcurl")
options(repos = BiocManager::repositories())

cat("R version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
print(BiocManager::repositories())

if (as.character(BiocManager::version()) != "3.18") {
  BiocManager::install(version = "3.18", ask = FALSE, lib = .libPaths()[1])
}

options(repos = BiocManager::repositories())

# ----------------------------
# CRAN packages
# ----------------------------
cran_pkgs <- c(
  "remotes",
  "fs",
  "s2",
  "dotCall64",
  "caTools",
  "patchwork",
  "dplyr",
  "ggplot2",
  "FNN",
  "slam",
  "deldir",
  "xgboost",
  "stringr",
  "RColorBrewer"
)

install_if_missing_cran(cran_pkgs)

# ----------------------------
# Core Bioconductor packages
# ----------------------------
core_bioc <- c(
  "limma",
  "BiocParallel",
  "sparseMatrixStats",
  "SingleCellExperiment",
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "chromVAR",
  "motifmatchr",
  "Gviz",
  "monocle3"
)

install_if_missing_bioc(core_bioc, force = FALSE, dependencies = TRUE)

# ----------------------------
# Genome packages
# ----------------------------
genome_pkgs <- c(
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10"
)

install_if_missing_bioc(genome_pkgs, force = FALSE, dependencies = TRUE)

# ----------------------------
# Heavy CRAN packages with explicit prereq checks
# ----------------------------
if (!is_installed("dotCall64")) {
  stop("dotCall64 is required before installing Seurat")
}
if (!is_installed("sparseMatrixStats")) {
  stop("sparseMatrixStats is required before installing Signac")
}

install_if_missing_cran("Seurat")
install_if_missing_cran("Signac")

# ----------------------------
# Diagnostics
# ----------------------------
cat("limma installed in primary lib:", is_installed("limma"), "\n")
cat("monocle3 installed in primary lib:", is_installed("monocle3"), "\n")
cat("chromVAR installed in primary lib:", is_installed("chromVAR"), "\n")
cat("motifmatchr installed in primary lib:", is_installed("motifmatchr"), "\n")

cat("limma path:\n")
print(tryCatch(find.package("limma"), error = function(e) NA_character_))

cat("monocle3 path:\n")
print(tryCatch(find.package("monocle3"), error = function(e) NA_character_))

# ----------------------------
# GitHub packages
# ----------------------------
install_github_if_missing(
  pkg_name = "cicero",
  repo = "cole-trapnell-lab/cicero-release",
  ref = "monocle3",
  upgrade = "never",
  dependencies = FALSE,
  force = FALSE
)

install_github_if_missing(
  pkg_name = "DIRECTNET",
  repo = "zhanglhbioinfor/DIRECT-NET",
  upgrade = "never",
  dependencies = FALSE,
  force = FALSE
)

cat("\nAll requested packages are installed and loadable.\n")