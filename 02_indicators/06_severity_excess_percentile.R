# ============================================================
# Excess and percentile-deviation severity
# ------------------------------------------------------------
#   HWES    sum of (Tmax - threshold) over heatwave days, in degree days
#   HWPD    sum of (estimated percentile - base) over heatwave days
#   EXCESS  sum of (precipitation - threshold) over extreme days, in mm
#   RPCTD   sum of (estimated percentile - base) over extreme days
#
# The percentile of each day is estimated by linear interpolation between the
# two percentile thresholds, capped at 100.
#
# Requires: steps 03 and 05
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))

yr        <- TARGET_YEAR
GDAL_OPTS <- c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6")

thr_TX90_file <- file.path(DIR_INDICES, "TX90p/Thresholds/TX90p_thresholds.tif")
thr_TX95_file <- file.path(DIR_INDICES, "TX95p/Thresholds/TX95p_thresholds.tif")
thr_R95_file  <- file.path(DIR_INDICES, "R95p/Thresholds/R95p_thresholds_wetday.tif")
thr_R99_file  <- file.path(DIR_INDICES, "R99p/Thresholds/R99p_thresholds_wetday.tif")

# PART 1: HEATWAVES (HWES + HWPD) 
cat("> Loading Tmax", yr, "+ thresholds ...\n")
tmax <- load_year("Tmax", yr); nD <- nlyr(tmax)        # daily Tmax in C
thr_TX90 <- rast(thr_TX90_file); thr_TX95 <- rast(thr_TX95_file)
pctile_slope <- clamp(5 / (thr_TX95 - thr_TX90), lower = 0, upper = 50)

make_hwes_fun <- function(nD) function(v) {
  exc <- v[1:nD]; t <- v[(nD+1):(2*nD)]; th <- v[(2*nD+1):(3*nD)]
  ok <- !is.na(exc) & (exc == 1); if (!any(ok)) return(0)
  r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  valid <- which(r$values & r$lengths >= 3); if (!length(valid)) return(0)
  idx <- unlist(Map(seq, starts[valid], ends[valid]), use.names = FALSE)
  sum(t[idx] - th[idx], na.rm = TRUE)
}
make_hwpd_fun <- function(nD, base) function(v) {
  exc <- v[1:nD]; t <- v[(nD+1):(2*nD)]; t90 <- v[(2*nD+1):(3*nD)]; sl <- v[(3*nD+1):(4*nD)]
  ok <- !is.na(exc) & (exc == 1); if (!any(ok)) return(0)
  r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  valid <- which(r$values & r$lengths >= 3); if (!length(valid)) return(0)
  idx <- unlist(Map(seq, starts[valid], ends[valid]), use.names = FALSE)
  p <- pmin(90 + sl[idx] * (t[idx] - t90[idx]), 100)
  sum(p - base, na.rm = TRUE)
}

hw_cfg <- list(
  list(tag = "TX90p", thr_file = thr_TX90_file, base = 90),
  list(tag = "TX95p", thr_file = thr_TX95_file, base = 95)
)
for (cfg in hw_cfg) {
  cat("> Heatwave severity", cfg$tag, "\n")
  thr <- rast(cfg$thr_file); exc <- tmax > thr[[1:nD]]
  res_hwes <- app(c(exc, tmax, thr[[1:nD]]), make_hwes_fun(nD))
  res_hwpd <- app(c(exc, tmax, thr_TX90[[1:nD]], pctile_slope[[1:nD]]),
                  make_hwpd_fun(nD, cfg$base))
  for (mn in c("HWES","HWPD")) {
    d <- file.path(DIR_INDICES, "Heatwaves", cfg$tag, "ANN", mn)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    r_out <- if (mn == "HWES") res_hwes else res_hwpd
    writeRaster(r_out, file.path(d, sprintf("%s_ANN_%s_%d.tif", cfg$tag, mn, yr)),
                datatype = "FLT4S", gdal = GDAL_OPTS, overwrite = TRUE)
    cat("  OK", mn, "\n")
  }
  rm(thr, exc, res_hwes, res_hwpd); gc()
}
rm(tmax, thr_TX90, thr_TX95, pctile_slope); gc()

# PART 2: RAINFALL (EXCESS + RPCTD) 
cat("> Loading precip", yr, "+ thresholds ...\n")
prec <- load_year("Precip", yr); nP <- nlyr(prec)       # daily precip in mm
thr_R95 <- rast(thr_R95_file); thr_R99 <- rast(thr_R99_file)
rain_slope <- clamp(4 / (thr_R99 - thr_R95), lower = 0, upper = 50)

make_excess_fun <- function(nP) function(v) {
  p <- v[1:nP]; b <- v[nP+1]; if (is.na(b)) return(0)
  ex <- which(!is.na(p) & p > b); if (!length(ex)) return(0)
  sum(p[ex] - b, na.rm = TRUE)
}
make_rpctd_fun <- function(nP, base) function(v) {
  p <- v[1:nP]; r95 <- v[nP+1]; sl <- v[nP+2]; b <- v[nP+3]
  if (is.na(b) || is.na(r95) || is.na(sl)) return(0)
  ex <- which(!is.na(p) & p > b); if (!length(ex)) return(0)
  pc <- pmin(95 + sl * (p[ex] - r95), 100)
  sum(pc - base, na.rm = TRUE)
}

rain_cfg <- list(
  list(tag = "R95p", thr_file = thr_R95_file, root = file.path(DIR_INDICES, "R95p"), base = 95),
  list(tag = "R99p", thr_file = thr_R99_file, root = file.path(DIR_INDICES, "R99p"), base = 99)
)
for (cfg in rain_cfg) {
  cat("> Rainfall severity", cfg$tag, "\n")
  base_thr <- rast(cfg$thr_file)
  res_excess <- app(c(prec, base_thr), make_excess_fun(nP))
  res_rpctd  <- app(c(prec, thr_R95, rain_slope, base_thr), make_rpctd_fun(nP, cfg$base))
  d_ex <- file.path(cfg$root, "ANNUAL_EXCESS");  dir.create(d_ex, recursive = TRUE, showWarnings = FALSE)
  d_pd <- file.path(cfg$root, "ANNUAL_PCTLDEV"); dir.create(d_pd, recursive = TRUE, showWarnings = FALSE)
  writeRaster(res_excess,
              file.path(d_ex, sprintf("%sEXCESS_%d.tif", cfg$tag, yr)),
              datatype = "FLT4S", gdal = GDAL_OPTS, overwrite = TRUE)
  writeRaster(res_rpctd,
              file.path(d_pd, sprintf("%sPCTD_%d.tif", cfg$tag, yr)),
              datatype = "FLT4S", gdal = GDAL_OPTS, overwrite = TRUE)
  cat("  Excess and percentile deviation written\n")
  rm(base_thr, res_excess, res_rpctd); gc()
}

message("Done. Next: step 07.")
