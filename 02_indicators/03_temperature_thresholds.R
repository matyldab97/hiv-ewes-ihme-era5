# ============================================================
# Day-of-year temperature thresholds
# ------------------------------------------------------------
# For each calendar day, the 90th and 95th percentile of daily maximum
# temperature across the reference years. The resulting 366-day series is
# smoothed with a centred 31-day circular moving average.
#
# Input:  data/era5/tmax
# Output: results/indices/TX90p/Thresholds, results/indices/TX95p/Thresholds
# Used by: steps 04 and 06
# ============================================================

source(file.path("01_setup", "01_config.R"))

suppressPackageStartupMessages(library(zoo))

# percentile sets to run (TX only)
PCTS <- list(`90p` = 0.90, `95p` = 0.95)

# ---- helpers ------------------------------------------------------------
smooth31 <- function(x) {                        # circular 31-day smoother
  padded <- c(tail(x, 15), x, head(x, 15))
  s <- rollapply(padded, 31, mean, align = "center", fill = NA, na.rm = TRUE)
  s[16:(15 + length(x))]
}
make_dirs <- function(root) {
  v <- file.path(root, "Thresholds")
  names(v) <- "thr"
  lapply(v, dir.create, recursive = TRUE, showWarnings = FALSE)
  as.list(v)
}

# ---- process each percentile ------------------------------------------------------------
for (pct_label in names(PCTS)) {
  pct  <- PCTS[[pct_label]]
  root <- file.path(DIR_INDICES, sprintf("TX%s", pct_label))
  dirs <- make_dirs(root)
  cat("\n=====  TX", pct_label, " =====\n", sep = "")

  # ---- threshold (cached) ------------------------------------------------------------
  thr_file <- file.path(dirs$thr, sprintf("TX%s_thresholds.tif", pct_label))
  if (file.exists(thr_file)) {
    cat("threshold cached - loading\n"); thr <- rast(thr_file)
  } else {
    cat("building 1991-2020 smoothed ", pct_label, " threshold\n")
    ref_all <- do.call(c, lapply(REF_YEARS, load_year, var = "Tmax"))
    doy <- unlist(lapply(REF_YEARS, \(y) if (is_leap_year(y)) 1:366 else 1:365))
    raw <- tapp(ref_all, doy, \(v) quantile(v, pct, na.rm = TRUE))
    thr <- app(raw, smooth31)
    writeRaster(thr, thr_file, datatype = "FLT4S", gdal = "COMPRESS=LZW", overwrite = TRUE)
    rm(ref_all); gc()
  }

}

message("Done. Next: step 04.")
