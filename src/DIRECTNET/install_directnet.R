args <- commandArgs(trailingOnly = TRUE)

personal_lib <- args[1]

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
  lib = personal_lib,
  upgrade = "never"
)
