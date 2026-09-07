# ============================================================
# Shared configuration
# ------------------------------------------------------------
# Sourced by every analysis script. Defines the ERA5 variable table, the
# reference period, and the daily loader that applies unit conversion.
#
# ERA5 file layout under DIR_ERA5:
#   pr/tmax/pet/ro/swvl*  yearly  files  era5_<token>_<YYYY>.nc
#
# ============================================================

source(file.path("R", "01_setup", "00_paths.R"))
suppressPackageStartupMessages({ library(terra) })

# ---- input variables ---------------------------------------------------
# convert() turns the raw ERA5 units into the units the indicators expect.
# A variable with no convert field is used as stored. Each is stored as one
# file per year, holding one layer per day.
ERA5_VARS <- list(
  Precip = list(dir = file.path(DIR_ERA5, "pr"),   token = "pr",
                convert = function(r) r * 1000),      # m -> mm
  Tmax   = list(dir = file.path(DIR_ERA5, "tmax"), token = "tmax",
                convert = function(r) r - 273.15),    # K -> degrees C
  PET    = list(dir = file.path(DIR_ERA5, "pet"),  token = "pet",
                convert = function(r) r * -1000),     # m -> mm; ERA5 stores potential
                                                      # evaporation as a negative
                                                      # downward flux, so the sign flips
  Runoff = list(dir = file.path(DIR_ERA5, "ro"),   token = "ro",
                convert = function(r) r * 1000),      # m -> mm
  Soil   = list(dir = file.path(DIR_RESULTS, "soil"), token = "soilwavg")
                                                      # m3/m3, dimensionless. Derived by
                                                      # step 02, not downloaded, so it is
                                                      # written under results rather than
                                                      # alongside the ERA5 inputs
)

# Raw soil layers and their thickness weights over the top 1 m.
# Read directly by step 02, which writes the depth-weighted mean as Soil.
SOIL_LAYERS <- list(
  L1 = list(dir = file.path(DIR_ERA5, "swvl1"), token = "swvl1", weight = 0.07),  #   0-7 cm
  L2 = list(dir = file.path(DIR_ERA5, "swvl2"), token = "swvl2", weight = 0.21),  #  7-28 cm
  L3 = list(dir = file.path(DIR_ERA5, "swvl3"), token = "swvl3", weight = 0.72)   # 28-100 cm
)

# ---- analysis window ---------------------------------------------------
REF_YEARS   <- 1991:2020
TARGET_YEAR <- 2018
CLAMP_BOUND <- 4     # bound on standardised drought indices, in SD

# The span the monthly drought series is built over. The distribution is always
# fitted on REF_YEARS; the index itself is computed across this whole span, and
# the twelve target-year layers are then taken out of it. The two coincide while
# the target year sits inside the reference period, which is the case here.
#
ANALYSIS_YEARS <- seq(min(REF_YEARS, TARGET_YEAR), max(REF_YEARS, TARGET_YEAR))

if (!(TARGET_YEAR %in% REF_YEARS))
  message(sprintf(paste("TARGET_YEAR (%d) lies outside REF_YEARS (%d-%d).",
                        "Fitting on REF_YEARS, series built over %d-%d;",
                        "every year in that span needs its input files."),
                  TARGET_YEAR, min(REF_YEARS), max(REF_YEARS),
                  min(ANALYSIS_YEARS), max(ANALYSIS_YEARS)))

# Values reaching the lower bound are separated into genuine extremes and cases
# where the baseline cannot support the fit at all because all values are close to 0. 
#
# Runoff degenerates by grid cell: over much of the domain the baseline is zero
# in most months all year round. 
SRI_MIN_BASELINE_SD <- 1.0   # mm, on the monthly baseline series
SRI_MAX_RATIO       <- 0.7   # target-year mean relative to its own baseline
#
# Precipitation degenerates by season: in an arid dry season the
# accumulation window is close to zero in every reference year.The ratio condition 
# that the runoff rule uses was tested for precipitation on the full domain and is
# not needed: at 0.7 it removes no cell-months at all. Changing the threshold 
# from 1 mm to 50 mm moves the reclassified population from 1.69%
# to 1.89% of PLHIV in the analysis domain, and every value from 10 mm upwards
# gives an identical set.
SPI_MIN_WINDOW_MM   <- 5     # mm, mean accumulation of the window over REF_YEARS

# ---- helpers -----------------------------------------------------------

is_leap_year <- function(y) (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)

# One variable-year as a daily SpatRaster (365 or 366 layers), converted.
era5_year_files <- function(v, year)
  file.path(v$dir, sprintf("era5_%s_%d.nc", v$token, year))

load_year <- function(var, year) {
  v <- ERA5_VARS[[var]]
  if (is.null(v))
    stop("Variable is not in ERA5_VARS, add it to R/01_setup/01_config.R: ", var)

  f <- era5_year_files(v, year)
  if (!file.exists(f)) stop("Missing ERA5 file: ", f)
  r <- rast(f)

  if (!is.null(v$convert)) r <- v$convert(r)

  exp_days <- if (is_leap_year(year)) 366 else 365
  if (nlyr(r) != exp_days)
    warning(sprintf("%s %d: %d layers (expected %d)", var, year, nlyr(r), exp_days))
  r
}

# ---- terra ------------------------------------------------------------
terraOptions(memfrac = as.numeric(Sys.getenv("TERRA_MEMFRAC", "0.4")),
             tempdir = DIR_TMP)

message("Configuration loaded. ERA5: ", DIR_ERA5, " | indicators: ", DIR_INDICES)
