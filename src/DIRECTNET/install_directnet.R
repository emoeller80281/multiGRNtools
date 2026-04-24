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
