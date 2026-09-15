## =============================================================================
## 00_inputs/01_era5_inputs.R
##
## Specifies the ERA5 input for the pipeline.The data are obtained from the
## Copernicus Climate Data Store and prepared once, outside this repository.
##
## INPUT   daily ERA5 NetCDF under data/era5/<var>/, prepared as described below
## OUTPUT  a report on the console
##
## -----------------------------------------------------------------------------
## WHERE THE DATA COME FROM
##
## Dataset     ERA5 hourly data on single levels from 1940 to present
##             (CDS identifier: reanalysis-era5-single-levels). Hersbach et al.
##             2020, Q J R Meteorol Soc 146:1999, doi 10.1002/qj.3803.
## Access      https://cds.climate.copernicus.eu, free account required.
##             Accepting the licence for the dataset is a prerequisite.
## Product     reanalysis
## Resolution  0.25 degrees, the native grid of the product
## Area        North 40, West -25, South -40, East 55, which covers Africa.
##             The CDS applies this subset server-side.
## Period      every year from 1991 to 2020 inclusive, all months, all days,
##             all 24 hours. The reference period and the target year are set in
##             01_setup/01_config.R; if either is changed, the years needed
##             here change with it, and load_year() names the first file it
##             cannot find.
##
## Variables, given by their CDS names:
##
##     2m_temperature                    daily maximum   -> data/era5/tmax
##     total_precipitation               daily sum       -> data/era5/pr
##     potential_evaporation             daily sum       -> data/era5/pet
##     runoff                            daily sum       -> data/era5/ro
##     volumetric_soil_water_layer_1     daily mean      -> data/era5/swvl1
##     volumetric_soil_water_layer_2     daily mean      -> data/era5/swvl2
##     volumetric_soil_water_layer_3     daily mean      -> data/era5/swvl3
##
## -----------------------------------------------------------------------------
## AGGREGATION FROM HOURLY TO DAILY
##
## The CDS produces hourly data. They are reduced to one value per day:
##
##     maximum   for temperature
##     sum       for precipitation, potential evaporation and runoff
##     mean      for the soil water layers, which are states rather than fluxes
##
## ERA5 stamps an accumulated field with the end of the hour it covers.The time
## coordinate of every accumulated variable is therefore shifted back by one
## hour, so that each hourly value is assigned to the day on which it actually
## fell. 
##
##
## -----------------------------------------------------------------------------
## LAYOUT AND FILE NAMES
##
## load_year() builds paths itself, so the names must match exactly:
##
##     data/era5/tmax/era5_tmax_YYYY.nc        one file per year
##     data/era5/pr/era5_pr_YYYY.nc            one file per year
##     data/era5/pet/era5_pet_YYYY.nc          one file per year
##     data/era5/ro/era5_ro_YYYY.nc            one file per year
##     data/era5/swvl1/era5_swvl1_YYYY.nc      one file per year, layers 1 to 3
##
##
## Each file holds one layer per day, 365 layers in a common year and 366 in a
## leap year, in calendar order.
##
##
## =============================================================================

source(file.path("01_setup", "01_config.R"))

suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0)

## Hour->day conversion is done outside R, so the check works on what reaches
## the pipeline: the files load_year() will open.
DOWNLOADED <- c("Precip", "Tmax", "PET", "Runoff")


## The same resolution load_year() uses, so the check reports exactly the files
## the pipeline will open.
expected_files <- function(var, year) era5_year_files(ERA5_VARS[[var]], year)

message("  ERA5 input check: ", min(ANALYSIS_YEARS), " to ", max(ANALYSIS_YEARS))

## ---- 1. every file the pipeline will ask for -------------------------------
missing <- character(0)
for (var in DOWNLOADED) {
  paths <- unlist(lapply(ANALYSIS_YEARS, function(y) expected_files(var, y)))
  gone  <- paths[!file.exists(paths)]
  message(sprintf("  %-7s %4d of %4d files present", var,
                  length(paths) - length(gone), length(paths)))
  missing <- c(missing, gone)
}

for (nm in names(SOIL_LAYERS)) {
  s <- SOIL_LAYERS[[nm]]
  paths <- unlist(lapply(ANALYSIS_YEARS, function(y) era5_year_files(s, y)))
  gone  <- paths[!file.exists(paths)]
  message(sprintf("  %-7s %4d of %4d files present", s$token,
                  length(paths) - length(gone), length(paths)))
  missing <- c(missing, gone)
}

if (length(missing)) {
  message("\n  Missing, first 10 of ", length(missing), ":")
  for (f in head(missing, 10)) message("    ", f)
  stop("ERA5 input is incomplete. Download and prepare the missing files as ",
       "described at the top of this script, then run it again.")
}

## ---- 2. layer counts and grid, on one year -------------------------------
## The observed value range of each variable is printed alongside, as a quick
## way to see that a file holds what its name says and is still in raw units.
message("\n  Checking layer counts and grid on ", TARGET_YEAR)
exp_days <- if (is_leap_year(TARGET_YEAR)) 366 else 365
problems <- character(0)
ref_geom <- NULL

check_one <- function(label, r) {
  if (nlyr(r) != exp_days)
    problems <<- c(problems, sprintf("%s: %d layers, expected %d",
                                     label, nlyr(r), exp_days))
  if (is.null(ref_geom)) {
    ref_geom <<- r[[1]]
  } else if (!compareGeom(r[[1]], ref_geom, stopOnError = FALSE)) {
    problems <<- c(problems, sprintf("%s: grid differs from the first variable read",
                                     label))
  }
  mm <- range(values(r[[seq_len(min(10, nlyr(r)))]], na.rm = TRUE), na.rm = TRUE)
  message(sprintf("    %-7s %3d layers | values %.4g to %.4g", label, nlyr(r), mm[1], mm[2]))
}

for (var in DOWNLOADED) {
  r <- rast(expected_files(var, TARGET_YEAR))
  check_one(var, r)
}

for (nm in names(SOIL_LAYERS)) {
  s <- SOIL_LAYERS[[nm]]
  r <- rast(era5_year_files(s, TARGET_YEAR))
  check_one(s$token, r)
}

if (length(problems)) {
  message("\n  Problems found:")
  for (p in problems) message("    ", p)
  stop("The ERA5 input does not have the expected shape.")
}

message("\n  Grid: ", paste(dim(ref_geom)[1:2], collapse = " x "), " cells at ",
        paste(round(res(ref_geom), 4), collapse = " by "), " degrees")
message("  ERA5 input complete and consistent.")
message("\nNext: 02_indicators/02_build_soil_moisture.R")
