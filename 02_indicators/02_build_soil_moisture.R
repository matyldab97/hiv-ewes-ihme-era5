# ============================================================
# Depth-weighted soil moisture
# ------------------------------------------------------------
# Combines ERA5 volumetric soil water layers 1 to 3 into a single daily series
# for the top 1 m, weighted by layer thickness (0.07, 0.21 and 0.72).
#
# Input:  data/era5/swvl1, swvl2, swvl3
# Output: results/soil
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))

out_dir <- ERA5_VARS$Soil$dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

w <- vapply(SOIL_LAYERS, `[[`, numeric(1), "weight")
w <- w / sum(w)
message(sprintf("Soil weights: L1=%.3f L2=%.3f L3=%.3f", w["L1"], w["L2"], w["L3"]))

# ANALYSIS_YEARS, so the target year is always built even when it falls outside
# the reference period. It equals REF_YEARS whenever the two coincide.
for (year in ANALYSIS_YEARS) {
  fL <- vapply(SOIL_LAYERS, function(L)
    file.path(L$dir, sprintf("era5_%s_%d.nc", L$token, year)),
    character(1))
  if (!all(file.exists(fL))) {
    message("skip ", year, " (missing layer file)"); next
  }

  r1 <- rast(fL["L1"]); r2 <- rast(fL["L2"]); r3 <- rast(fL["L3"])
  if (!(nlyr(r1) == nlyr(r2) && nlyr(r2) == nlyr(r3)))
    stop("Layer day-count mismatch at ", year)

  soil <- w["L1"] * r1 + w["L2"] * r2 + w["L3"] * r3   # daily weighted mean

  outf <- file.path(out_dir, sprintf("era5_%s_%d.nc", ERA5_VARS$Soil$token, year))
  writeCDF(soil, outf, overwrite = TRUE, varname = "soilwavg",
           longname = "depth-weighted soil moisture (top 1m)", unit = "m3 m-3")
  message("done ", year)
}

message("Done. Next: step 03.")
