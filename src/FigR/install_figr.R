args <- commandArgs(trailingOnly = TRUE)

personal_lib <- args[1]

remotes::install_github(
  "buenrostrolab/FigR",
  lib = personal_lib,
  upgrade = "never"
)