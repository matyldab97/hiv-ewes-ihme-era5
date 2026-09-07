# ============================================================
# Standardised Precipitation Index
# ------------------------------------------------------------
# Monthly precipitation totals fitted to a gamma distribution over the
# reference period, at 1 and 3 month accumulation, after McKee, Doesken and
# Kleist (1993, Proceedings of the 8th Conference on Applied Climatology, 179).
#
#
# Because the monthly index is what step 11 reads, the severity sums and the
# month counts both follow from this one operation.
#
# Output: results/indices/SPI
# ============================================================

source(file.path("R", "01_setup", "01_config.R")); source(file.path("R", "01_setup", "02_lib_indices.R"))

out_dir <- file.path(DIR_INDICES, "SPI"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
SCALES <- c(1, 3)
THRS    <- c(-1.0, -1.5)

monthly <- build_monthly_ref("Precip", "sum")          # 360 layers YYYY_MM

yrs <- as.integer(substr(names(monthly), 1, 4))
mos <- as.integer(substr(names(monthly), 6, 7))

# Raw accumulation over the scale-month window ending in each month, in mm.
window_accumulation <- function(scale) {
  if (scale == 1) return(monthly)
  roll <- app(monthly, function(v) zoo::rollapply(v, width = scale, FUN = sum,
                                                 na.rm = TRUE, fill = NA, align = "right"))
  names(roll) <- names(monthly)
  roll
}

# Mean of that accumulation over the reference years, one layer per calendar
# month, aligned with the twelve target-year index layers.
window_norm <- function(acc) {
  nrm <- rast(lapply(1:12, function(m)
    app(acc[[which(yrs %in% REF_YEARS & mos == m)]], mean, na.rm = TRUE)))
  names(nrm) <- sprintf("%d_%02d", TARGET_YEAR, 1:12)
  nrm
}

for (sc in SCALES) {
  message("SPI-", sc)
  idx  <- standardised_stack(monthly, sc, kind = "spi")
  yr   <- target_year_layers(idx)                       # 12 layers, target year

  nrm <- window_norm(window_accumulation(sc))

  clamped    <- !is.na(yr) & yr <= -CLAMP_BOUND + 1e-9
  unfittable <- clamped & (nrm < SPI_MIN_WINDOW_MM)
  message(sprintf(paste("    SPI-%d: %d cell-months at the lower bound =",
                        "%d unfittable (window norm below %g mm) + %d genuine"),
                  sc, sum(values(clamped), na.rm = TRUE),
                  sum(values(unfittable), na.rm = TRUE), SPI_MIN_WINDOW_MM,
                  sum(values(clamped), na.rm = TRUE) - sum(values(unfittable), na.rm = TRUE)))

  yr <- ifel(unfittable, 0, yr)
  names(yr) <- sprintf("%d_%02d", TARGET_YEAR, 1:12)

  if (sc == 3)
    writeRaster(app(unfittable, function(v) as.numeric(any(v > 0, na.rm = TRUE))),
                file.path(out_dir, sprintf("SPI3_%d_unfittable_cells.tif", TARGET_YEAR)),
                overwrite = TRUE)

  writeRaster(yr,
              file.path(out_dir, sprintf("SPI%d_%d_monthly.tif", sc, TARGET_YEAR)),
              overwrite = TRUE)
  for (thr in THRS) {
    writeRaster(severity_sum(yr, thr),
                file.path(out_dir, sprintf("SPI%d_%d_drought_severity_sum_leq%.1f.tif",
                                           sc, TARGET_YEAR, thr)), overwrite = TRUE)
  }
}

message("Done. Next: step 08.")
