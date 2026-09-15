# ============================================================
# Heatwave metrics
# ------------------------------------------------------------
# A heatwave is a run of at least three consecutive days on which daily maximum
# temperature exceeds the day-of-year threshold from step 03.
#
#   HWN   number of events
#   HWF   total heatwave days
#   HWM   mean daily maximum temperature across heatwave days
#   HWMF  HWM x HWF, the cumulative heatwave temperature
#   HWNM  HWM x HWN
#
# HWMF is the metric reported. The other four are written because step 13
# screens all of them side by side, and because HWF is the day count the
# threshold exposure in step 15 is defined on.
#
# The family follows Perkins and Alexander (2013, Journal of Climate 26:4500,
# doi 10.1175/JCLI-D-12-00383.1). A three-day minimum is used rather than the
# six days of the ETCCDI definition, because shorter events already carry
# health effects and the shorter minimum is the one used in vulnerability
# frameworks for this region.
#
# Output: results/indices/Heatwaves/TX90p, results/indices/Heatwaves/TX95p
# ============================================================

source(file.path("01_setup", "01_config.R"))

TAGS <- list(
  TX90p = file.path(DIR_INDICES, "TX90p/Thresholds/TX90p_thresholds.tif"),
  TX95p = file.path(DIR_INDICES, "TX95p/Thresholds/TX95p_thresholds.tif")
)
out_root <- file.path(DIR_INDICES, "Heatwaves")

hw_stats <- function(is_hw_day, temp) {
  ok <- !is.na(is_hw_day) & (is_hw_day == 1)
  if (!any(ok)) return(c(0, 0, NA_real_, 0, 0))
  r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  valid <- which(r$values & r$lengths >= 3)
  HWN <- length(valid); if (HWN == 0) return(c(0, 0, NA_real_, 0, 0))
  HWF <- sum(r$lengths[valid])
  idx <- unlist(Map(seq, starts[valid], ends[valid]), use.names = FALSE)
  HWM <- mean(temp[idx], na.rm = TRUE)
  c(HWN, HWF, HWM, HWM * HWF, HWM * HWN)
}

tday0 <- load_year("Tmax", TARGET_YEAR); nD <- nlyr(tday0)
metrics <- c("HWN", "HWF", "HWM", "HWMF", "HWNM")

for (tag in names(TAGS)) {
  cat("->", tag, "\n")
  thr  <- rast(TAGS[[tag]])
  days <- if (tag == names(TAGS)[1]) tday0 else load_year("Tmax", TARGET_YEAR)
  exc  <- days > thr[[1:nD]]
  st   <- c(exc, days)
  res  <- app(st, fun = function(v) hw_stats(v[1:nD] == 1, v[(nD+1):(2*nD)]))
  names(res) <- metrics
  for (mn in metrics) {
    d <- file.path(out_root, tag, "ANN", mn); dir.create(d, showWarnings = FALSE, recursive = TRUE)
    writeRaster(res[[mn]],
                file.path(d, sprintf("%s_ANN_%s_%d.tif", tag, mn, TARGET_YEAR)),
                overwrite = TRUE)
  }
}

message("Done. Next: step 05.")
