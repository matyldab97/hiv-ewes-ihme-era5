# ============================================================
# Path configuration
# ------------------------------------------------------------
# Every path in this repository resolves from a single root. Set it once,
# either as an environment variable
#     Sys.setenv(EWE_ROOT = "/path/to/project")
# or in a local, untracked file R/01_setup/00_paths_local.R that assigns
# EWE_ROOT. Nothing else in the repository contains an absolute path.
#
# Expected layout under EWE_ROOT:
#
#   data/era5/<var>/            daily ERA5 NetCDF, one folder per variable
#   data/plhiv/                 IHME gridded PLHIV counts, 5x5km grids
#   data/prevalence/            IHME gridded HIV prevalence, 5x5km grids
#   data/boundaries/gadm.gpkg   GADM polygons, subset to the 43 SSA countries
#   results/soil/               depth-weighted soil moisture, built by step 02
#   results/indices/            extreme weather indicator rasters
#   results/montecarlo/         PLHIV draws and uncertainty intervals
#   results/exposure/           exposure classification tables and figures
#   results/screening/          indicator screening panels
#   results/gistar/inputs/      Gi* input stacks
#   results/gistar/outputs/     Gi* rasters and tables
#   results/gistar/plots/       Gi* maps
#   results/supplementary/      the day threshold sweep, step 20
#   tmp/                        scratch space for terra
#
# Sourced by 01_config.R. Not run directly.
# ============================================================

if (file.exists(file.path("R", "01_setup", "00_paths_local.R"))) {
  source(file.path("R", "01_setup", "00_paths_local.R"))
}

if (!exists("EWE_ROOT")) {
  EWE_ROOT <- Sys.getenv("EWE_ROOT", unset = NA_character_)
}
if (is.na(EWE_ROOT) || !nzchar(EWE_ROOT)) {
  stop("EWE_ROOT is not set. Either Sys.setenv(EWE_ROOT = \"/path/to/project\") ",
       "or create R/01_setup/00_paths_local.R containing ",
       "EWE_ROOT <- \"/path/to/project\".")
}
EWE_ROOT <- normalizePath(EWE_ROOT, mustWork = FALSE)

# An override is used when it is set and non-empty, whether it came from the
# environment or from 00_paths_local.R; otherwise the default under EWE_ROOT
# applies.
from_env <- function(name, default) {
  v <- if (exists(name, inherits = TRUE)) get(name, inherits = TRUE) else Sys.getenv(name)
  if (is.null(v) || is.na(v) || !nzchar(v)) default else normalizePath(v, mustWork = FALSE)
}

# ---- inputs ------------------------------------------------------------
DIR_DATA       <- from_env("EWE_DATA", file.path(EWE_ROOT, "data"))
DIR_ERA5       <- file.path(DIR_DATA, "era5")
DIR_PLHIV_FINE <- from_env("EWE_PLHIV", file.path(DIR_DATA, "plhiv"))
DIR_PREV_FINE  <- from_env("EWE_PREV",  file.path(DIR_DATA, "prevalence"))
FILE_GADM      <- from_env("EWE_GADM",  file.path(DIR_DATA, "boundaries", "gadm.gpkg"))

# ---- outputs -----------------------------------------------------------
DIR_RESULTS    <- from_env("EWE_RESULTS", file.path(EWE_ROOT, "results"))
DIR_INDICES    <- file.path(DIR_RESULTS, "indices")
DIR_MC         <- file.path(DIR_RESULTS, "montecarlo")
DIR_EXPOSURE   <- file.path(DIR_RESULTS, "exposure")
DIR_EXP_FIGS   <- file.path(DIR_EXPOSURE, "figures")  
DIR_EXP_TABLES <- file.path(DIR_EXPOSURE, "tables")   
DIR_SCREENING  <- file.path(DIR_RESULTS, "screening")
DIR_GI_IN      <- file.path(DIR_RESULTS, "gistar", "inputs")
DIR_GI_OUT     <- file.path(DIR_RESULTS, "gistar", "outputs")
DIR_GI_PLOTS   <- file.path(DIR_RESULTS, "gistar", "plots")
DIR_SUPP       <- file.path(DIR_RESULTS, "supplementary")
DIR_TMP        <- from_env("EWE_TMP", file.path(EWE_ROOT, "tmp"))

for (.d in c(DIR_RESULTS, DIR_INDICES, DIR_MC, DIR_EXPOSURE,
             DIR_EXP_FIGS, DIR_EXP_TABLES, DIR_SCREENING,
             DIR_GI_IN, DIR_GI_OUT, DIR_GI_PLOTS, DIR_SUPP, DIR_TMP)) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d, from_env)
