# ============================================================
# Drought month counts
# ------------------------------------------------------------
# Annual number of months at or below each drought threshold, for every index
# and accumulation scale.
#
# Requires: steps 07 to 10
# ============================================================

source(file.path("01_setup", "01_config.R"))

yr     <- TARGET_YEAR
SCALES <- c(1, 3)
THRS   <- c(-1.0, -1.5)
INDEX_DIRS <- list(
  SPI   = file.path(DIR_INDICES, "SPI"),
  SPEI  = file.path(DIR_INDICES, "SPEI"),
  SRI   = file.path(DIR_INDICES, "SRI"),
  SMA = file.path(DIR_INDICES, "SMA")
)

for (nm in names(INDEX_DIRS)) {
  d <- INDEX_DIRS[[nm]]
  for (sc in SCALES) {
    f <- file.path(d, sprintf("%s%d_%d_monthly.tif", nm, sc, yr))
    if (!file.exists(f)) { message("skip (missing): ", f); next }
    stk <- rast(f)                                  # 12 monthly index layers
    for (thr in THRS) {
      cnt <- app(stk <= thr, sum, na.rm = TRUE)     # months in drought
      outf <- file.path(d, sprintf("%s%d_%d_drought_month_count_leq%.1f.tif",
                                   nm, sc, yr, thr))
      writeRaster(cnt, outf, overwrite = TRUE)
      message(nm, sc, "  count<=", thr, "  -> ", basename(outf))
    }
  }
}

message("Done. Next: step 12.")
