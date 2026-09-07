# ============================================================
# Wet-day extreme rainfall indices
# ------------------------------------------------------------
# A wet day has at least 1 mm of precipitation. The threshold is the 95th or
# 99th percentile of the wet-day distribution over the reference period.
#
#   DAY  number of extreme rainfall days in the target year
#   TOT  precipitation accumulated on those days
#
# Output: results/indices/R95p, results/indices/R99p
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))

PCTS <- c("95" = 0.95, "99" = 0.99)

mk <- function(tag) {
  base <- file.path(DIR_INDICES, sprintf("R%sp", tag))
  p <- file.path(base, c("Thresholds","ANNUAL_DAY","ANNUAL_TOT"))
  names(p) <- c("thr","ann_day","ann_tot")
  lapply(p, dir.create, recursive = TRUE, showWarnings = FALSE); p
}

# reference daily stack (1991-2020) for thresholds
message("Reading 1991-2020 daily precip for thresholds ...")
ref_all <- do.call(c, lapply(REF_YEARS, load_year, var = "Precip"))

for (tag in names(PCTS)) {
  p <- mk(tag); pr <- PCTS[[tag]]
  thr_file <- file.path(p["thr"], sprintf("R%sp_thresholds_wetday.tif", tag))
  if (file.exists(thr_file)) { thr <- rast(thr_file) } else {
    message("Threshold R", tag, "p (wet-day ", pr*100, "th pctile)")
    thr <- app(ref_all, function(v) {
      w <- v[v >= 1 & !is.na(v)]; if (!length(w)) return(NA_real_)
      quantile(w, pr, na.rm = TRUE, names = FALSE)
    })
    writeRaster(thr, thr_file, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
  }

  yr <- TARGET_YEAR; rr <- load_year("Precip", yr)

  exc <- rr > thr                 # extreme day = above wet-day pctile
  amt <- rr; amt[!exc] <- 0                     # precip on extreme days

  cnt_ann <- app(exc, sum, na.rm = TRUE)
  tot_ann <- app(amt, sum, na.rm = TRUE)
  writeRaster(cnt_ann, file.path(p["ann_day"], sprintf("R%spDAY_%d.tif", tag, yr)), overwrite=TRUE)
  writeRaster(tot_ann, file.path(p["ann_tot"], sprintf("R%spTOT_%d.tif", tag, yr)), overwrite=TRUE)

  cat("R", tag, "p done.\n", sep = "")
}

message("Done. Next: step 06.")
