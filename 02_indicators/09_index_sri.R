# ============================================================
# Standardised Runoff Index
# ------------------------------------------------------------
# Monthly runoff totals fitted to a gamma distribution over the reference
# period.
#
#
# Output: results/indices/SRI
# ============================================================

source(file.path("01_setup", "01_config.R")); source(file.path("01_setup", "02_lib_indices.R"))

SRI_YEARS <- ANALYSIS_YEARS   # monthly series span; the distribution is always
                              # fitted on REF_YEARS, so this may be widened
out_dir <- file.path(DIR_INDICES, "SRI"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
SCALES <- c(1, 3); THRS <- c(-1.0, -1.5)

monthly <- build_monthly_ref("Runoff", "sum", years = SRI_YEARS)

baseline_sd   <- app(monthly, function(v) sd(v, na.rm = TRUE))
baseline_mean <- app(monthly, function(v) mean(v, na.rm = TRUE))
target_mean   <- app(target_year_layers(monthly), function(v) mean(v, na.rm = TRUE))
ratio         <- target_mean / baseline_mean

for (sc in SCALES) {
  message("SRI-", sc)
  idx <- standardised_stack(monthly, sc, kind = "spi")
  yr  <- target_year_layers(idx)

  # Cells at the lower bound in at least one month of the target year.
  clamped <- app(yr, function(v) as.numeric(any(!is.na(v) & v <= -CLAMP_BOUND + 1e-9)))

  genuine <- (clamped == 1) &
             (baseline_sd >= SRI_MIN_BASELINE_SD) &
             (ratio < SRI_MAX_RATIO)
  unfittable <- (clamped == 1) & !genuine

  n_clamped    <- sum(values(clamped) == 1, na.rm = TRUE)
  n_genuine    <- sum(values(genuine), na.rm = TRUE)
  n_unfittable <- sum(values(unfittable), na.rm = TRUE)
  message(sprintf("    SRI-%d: %d cells at the lower bound = %d genuine + %d unfittable",
                  sc, n_clamped, n_genuine, n_unfittable))

  yr[unfittable] <- 0

  if (sc == 3)
    writeRaster(unfittable,
                file.path(out_dir, sprintf("SRI3_%d_unfittable_cells.tif", TARGET_YEAR)),
                overwrite = TRUE)

  writeRaster(yr,
              file.path(out_dir, sprintf("SRI%d_%d_monthly.tif", sc, TARGET_YEAR)),
              overwrite = TRUE)
  for (thr in THRS)
    writeRaster(severity_sum(yr, thr),
      file.path(out_dir, sprintf("SRI%d_%d_drought_severity_sum_leq%.1f.tif",
                                 sc, TARGET_YEAR, thr)), overwrite = TRUE)
}

message("Done. Next: step 10.")
