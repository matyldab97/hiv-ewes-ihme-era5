# ============================================================
# Standardised Precipitation-Evapotranspiration Index
# ------------------------------------------------------------
# The climatic water balance (precipitation minus potential evapotranspiration)
# fitted to a log-logistic distribution over the reference period, after
# Vicente-Serrano, Begueria and Lopez-Moreno (2010, Journal of Climate 23:1696,
# doi 10.1175/2009JCLI2909.1).
#
#
# Output: results/indices/SPEI
# ============================================================

source(file.path("01_setup", "01_config.R")); source(file.path("01_setup", "02_lib_indices.R"))

out_dir <- file.path(DIR_INDICES, "SPEI"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
SCALES <- c(1, 3); THRS <- c(-1.0, -1.5)

prec <- build_monthly_ref("Precip", "sum")
pet  <- build_monthly_ref("PET",    "sum")
bal  <- prec - pet; names(bal) <- names(prec)        # monthly water balance

for (sc in SCALES) {
  message("SPEI-", sc)
  idx <- standardised_stack(bal, sc, kind = "spei")
  yr  <- target_year_layers(idx)
  writeRaster(yr,
              file.path(out_dir, sprintf("SPEI%d_%d_monthly.tif", sc, TARGET_YEAR)),
              overwrite = TRUE)
  for (thr in THRS)
    writeRaster(severity_sum(yr, thr),
      file.path(out_dir, sprintf("SPEI%d_%d_drought_severity_sum_leq%.1f.tif",
                                 sc, TARGET_YEAR, thr)), overwrite = TRUE)
}

message("Done. Next: step 09.")
