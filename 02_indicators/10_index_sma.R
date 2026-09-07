# ============================================================
# Soil Moisture Anomaly
# ------------------------------------------------------------
# For each grid cell and calendar month, the empirical distribution of the
# reference years is used to obtain the percentile rank of the observed value,
# which is then transformed to a standard normal deviate.
#
# Requires: step 02
# Output:   results/indices/SMA
# ============================================================

source(file.path("R", "01_setup", "01_config.R")); source(file.path("R", "01_setup", "02_lib_indices.R"))

# ---- SMA: empirical standardisation per calendar month ------------------
# Returns the 12 target-year monthly SMA layers for a given scale.
#
# Soil moisture is bounded and strongly seasonal, so no parametric family fits
# it across the year. Each calendar month is therefore ranked against the same
# calendar month of the reference years, ecdf then qnorm.
#
# The empirical rank cannot exceed the reference sample, so with 30 reference
# years the most extreme dry month reachable is qnorm(1/30) = -1.83. SMA
# therefore never reaches the -2 or -3 that SPI can, and its thresholds are not
# interchangeable with the parametric indices.
sma_target_year <- function(monthly_stack, scale,
                              ref_years = REF_YEARS, target_year = TARGET_YEAR) {
  roll <- app(monthly_stack, function(v)
    zoo::rollapply(v, width = scale, FUN = mean, na.rm = TRUE,
                   fill = NA, align = "right"))
  names(roll) <- names(monthly_stack)
  yrs <- as.integer(substr(names(roll), 1, 4))
  mos <- as.integer(substr(names(roll), 6, 7))
  out <- vector("list", 12)
  for (m in 1:12) {
    iref <- which(yrs %in% ref_years & mos == m)
    icur <- which(yrs == target_year   & mos == m)
    base_br <- roll[[iref]]; N <- nlyr(base_br)
    all_br  <- c(base_br, roll[[icur]])
    out[[m]] <- app(all_br, function(x) {
      if (all(is.na(x))) return(NA_real_)
      e <- ecdf(x[1:N]); qnorm(e(x[N + 1]))
    })
  }
  r <- rast(out); names(r) <- sprintf("%d_%02d", target_year, 1:12)
  app(r, clamp_uniform)            # uniform clamp, same as the other indices
}

SOIL_YEARS <- ANALYSIS_YEARS
out_dir <- file.path(DIR_INDICES, "SMA"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
SCALES <- c(1, 3); THRS <- c(-1.0, -1.5)

monthly <- build_monthly_ref("Soil", "mean", years = SOIL_YEARS)   # monthly MEAN

for (sc in SCALES) {
  message("SMA-", sc)
  yr <- sma_target_year(monthly, sc, ref_years = REF_YEARS,
                          target_year = TARGET_YEAR)        # 12 target-year layers
  writeRaster(yr,
              file.path(out_dir, sprintf("SMA%d_%d_monthly.tif", sc, TARGET_YEAR)),
              overwrite = TRUE)
  for (thr in THRS)
    writeRaster(severity_sum(yr, thr),
      file.path(out_dir, sprintf("SMA%d_%d_drought_severity_sum_leq%.1f.tif",
                                 sc, TARGET_YEAR, thr)), overwrite = TRUE)
}

message("Done. Next: step 11.")
