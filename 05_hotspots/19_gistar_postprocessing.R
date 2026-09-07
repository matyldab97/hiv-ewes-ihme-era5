# ============================================================
# Hotspot overlap
# ------------------------------------------------------------
# Cells where an HIV prevalence hotspot coincides with an indicator hotspot,
# and the number of PLHIV living in them.
#
# Requires: steps 12 and 18
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))
source(file.path("R", "01_setup", "03_lib_shared.R"))


# Reads the categorical and z-score rasters from step 18 and reports, both
# continent-wide and per country:
#
#   overlap maps            per indicator, cells that are hot in both the
#                           prevalence and the indicator surface, and cells
#                           that are cold in both
#   PLHIV counts            number and proportion of PLHIV in each category,
#                           at the 90, 95 and 99 per cent levels
#
# Main indicators: SPI3, TX90p, R95p. Six further indicators are reported as
# a sensitivity analysis.

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# ============================================================
# CONFIG
# ============================================================

config <- list(
  gi_dir = DIR_GI_OUT,
  cat_suffix = "_GiStar_categorical.tif",
  res_suffix = "_GiStar_results.tif",
  plot_dir = file.path(DIR_GI_PLOTS, "overlap_vulnerability"),

  # PLHIV source for population counts
  plhiv_ui_path = file.path(DIR_MC, sprintf("ALL_15_59_%d_UI_0p25.tif", TARGET_YEAR)),
  plhiv_draws_path = file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", TARGET_YEAR)),

  fdr_method = "BH",
  hiv_var    = "PREVpct_mean",

  gadm_path = FILE_GADM,

  main_vars = c(
    "TX90p_HWMF",
    "R95pTOT",
    "SPI3_leq1p0"
  ),

  robustness_vars = c(
    "TX95p_HWMF",
    "R99pTOT",
    "SPI3_leq1p5",
    "SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0"
  ),

  plot_width  = 14,
  plot_height = 8,

  var_labels = c(
    TX90p_HWMF    = "Heatwave (90th percentile)",
    R95pTOT       = "Extreme Rainfall (95th percentile)",
    SPI3_leq1p0   = "Drought (SPI, \u2264-1.0)",
    TX95p_HWMF    = "Heatwave (95th percentile)",
    R99pTOT       = "Extreme Rainfall (99th percentile)",
    SPI3_leq1p5   = "Drought (SPI, \u2264-1.5)",
    SPEI3_leq1p0  = "Drought (SPEI, \u2264-1.0)",
    SRI3_leq1p0   = "Drought (SRI, \u2264-1.0)",
    SMA3_leq1p0   = "Drought (SMA, \u2264-1.0)"
  )
)

for (d in c(config$plot_dir,
            file.path(config$plot_dir, "main"),
            file.path(config$plot_dir, "robustness")))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

message("  Gi* post-processing: hotspot overlap maps")

# ============================================================
# Load boundaries
# ============================================================

message("  Loading GADM boundaries...")
borders_v <- NULL
if (file.exists(config$gadm_path)) {
  borders_sf <- st_read(config$gadm_path, quiet = TRUE)
  borders_v <- vect(borders_sf)
  message(sprintf("    Loaded: %d country polygons", nrow(borders_v)))
}

# ============================================================
# Load PLHIV raster + 1000 draws
# ============================================================

plhiv_med <- NULL
plhiv_mat <- NULL
country_v <- NULL

if (file.exists(config$plhiv_ui_path)) {
  message("  Loading PLHIV UI raster...")
  ui_r <- rast(config$plhiv_ui_path)
  plhiv_idx <- grep("PLHIV_mean|PLHIV_median|PLHIV_med", names(ui_r))
  country_idx <- grep("^country$", names(ui_r))
  if (length(plhiv_idx) > 0) plhiv_med <- values(ui_r[[plhiv_idx[1]]], mat = FALSE)
  if (length(country_idx) > 0) country_v <- values(ui_r[[country_idx[1]]], mat = FALSE)
  rm(ui_r)
}

if (file.exists(config$plhiv_draws_path)) {
  message("  Loading 1000 PLHIV draws for CI...")
  draws_r <- rast(config$plhiv_draws_path)
  draw_idx <- grep("^PLHIV_d[0-9]", names(draws_r))
  if (length(draw_idx) == 0) draw_idx <- grep("^d[0-9]", names(draws_r))
  if (length(draw_idx) > 0) {
    plhiv_mat <- values(draws_r[[draw_idx]], mat = TRUE)
    plhiv_mat[!is.finite(plhiv_mat)] <- 0
    message(sprintf("    PLHIV draws: %d cells x %d draws", nrow(plhiv_mat), ncol(plhiv_mat)))
    if (is.null(plhiv_med)) plhiv_med <- rowMeans(plhiv_mat, na.rm = TRUE)
    # Get country from draws file if not from UI
    if (is.null(country_v)) {
      ci <- grep("^country$", names(draws_r))
      if (length(ci) > 0) country_v <- values(draws_r[[ci[1]]], mat = FALSE)
    }
  }
  rm(draws_r); gc(FALSE)
} else {
  message("  Warning: PLHIV draws file not found, reporting the mean without intervals")
}

if (!is.null(plhiv_med)) {
  message(sprintf("    PLHIV total: %s", format(sum(plhiv_med, na.rm = TRUE), big.mark = ",")))
}

draw_borders <- function(mode = "SSA") {
  if (!is.null(borders_v)) {
    lwd <- if (mode == "COUNTRY") 0.8 else 0.5
    lines(borders_v, col = "black", lwd = lwd)
  }
}

# ============================================================
# FIND FILES
# ============================================================

cat_esc <- gsub("\\.", "\\\\.", config$cat_suffix)
res_esc <- gsub("\\.", "\\\\.", config$res_suffix)

cat_files <- list.files(config$gi_dir, pattern = paste0(cat_esc, "$"), full.names = TRUE)
res_files <- list.files(config$gi_dir, pattern = paste0(res_esc, "$"), full.names = TRUE)

cat_ssa <- grep("_SSA", cat_files, value = TRUE)
cat_cty <- grep("_COUNTRY", cat_files, value = TRUE)
res_ssa <- grep("_SSA", res_files, value = TRUE)
res_cty <- grep("_COUNTRY", res_files, value = TRUE)

# ============================================================
# Process one mode
# ============================================================

process_mode <- function(cat_fp, res_fp, ref) {
  message(sprintf("\n  === %s MODE ===", ref))

  r_cat <- rast(cat_fp)
  r_res <- rast(res_fp)
  id <- sub(paste0(cat_esc, "$"), "", basename(cat_fp))
  fdr <- config$fdr_method

  # Everything this step reports is labelled as BH-corrected.
  
  cat_layer_name <- function(var) {
    nm <- paste0("CatFDR_", fdr, "_", var)
    if (nm %in% names(r_cat)) return(nm)
    if (paste0("Cat_", var) %in% names(r_cat))
      stop("The categorical raster carries only the uncorrected layer Cat_", var,
           ", not ", nm, ". This step reports ", fdr, "-corrected results, so the ",
           "uncorrected layer must not be used in its place. Rerun step 18, which ",
           "writes both.")
    NA_character_
  }
  has_cat_layer    <- function(var) !is.na(cat_layer_name(var))
  present_cat_vars <- function(vars) vars[vapply(vars, has_cat_layer, logical(1))]

  # HIV layer
  hiv_nm <- cat_layer_name(config$hiv_var)
  if (is.na(hiv_nm))
    stop("The categorical raster has no layer for the HIV surface (",
         config$hiv_var, "), so no overlap can be computed. Rerun step 18.")
  hiv_cat <- r_cat[[hiv_nm]]
  hiv_v <- values(hiv_cat, mat = FALSE)

  all_vars <- c(config$main_vars, config$robustness_vars)
  total_plhiv_ssa <- if (!is.null(plhiv_med)) sum(plhiv_med, na.rm = TRUE) else NULL

  # Threshold definitions: cat >= threshold_cat means hotspot at that level
  # Cat coding: 1 = 90%, 2 = 95%, 3 = 99%
  thresholds <- list(
    "90" = list(cat_min = 1, label = "90%"),
    "95" = list(cat_min = 2, label = "95%"),
    "99" = list(cat_min = 3, label = "99%")
  )
  main_threshold <- "95"  # default

  # Helper: make hotspot-only overlap raster at given threshold
  make_overlap <- function(ewe_v, hiv_v, thr_cat) {
    valid <- !is.na(ewe_v) & !is.na(hiv_v)
    ov <- rep(NA_integer_, length(ewe_v))
    ov[valid] <- 0L
    ov[valid & ewe_v >= thr_cat & hiv_v >= thr_cat] <- 1L  # both hotspot
    ov
  }

  col_ov <- c("#F5F0E1", "red")
  ov_brks <- c(-0.5, 0.5, 1.5)
  ov_leg <- c("Not an overlapping hotspot", "An overlapping hotspot")


  # ============================================================
  # 1. INDIVIDUAL OVERLAP MAPS + PLHIV (main threshold = 95%)
  # ============================================================

  overlap_rows <- list()

  for (var in all_vars) {
    ewe_nm <- cat_layer_name(var)
    if (!ewe_nm %in% names(r_cat)) { message(sprintf("    %s: skipped", var)); next }
    ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)

    thr_cat <- thresholds[[main_threshold]]$cat_min
    ov <- make_overlap(ewe_v, hiv_v, thr_cat)

    n_both <- sum(ov == 1, na.rm = TRUE)
    plhiv_both <- if (!is.null(plhiv_med)) round(sum(plhiv_med[ov == 1 & !is.na(ov)], na.rm = TRUE)) else NA
    pct_plhiv <- if (!is.na(plhiv_both) && !is.null(total_plhiv_ssa) && total_plhiv_ssa > 0) {
      round(100 * plhiv_both / total_plhiv_ssa, 2)
    } else NA

    is_main <- var %in% config$main_vars
    subdir <- if (is_main) "main" else "robustness"

    overlap_rows[[var]] <- data.frame(
      mode = ref, variable = var, label = get_label(var),
      analysis = if (is_main) "main" else "robustness",
      threshold = thresholds[[main_threshold]]$label,
      n_both_hotspot = n_both,
      plhiv_both_hotspot = plhiv_both, pct_plhiv_ssa = pct_plhiv,
      stringsAsFactors = FALSE)

    r_ov <- rast(hiv_cat); values(r_ov) <- ov

    jpeg(file.path(config$plot_dir, subdir,
                  paste0(id, "_overlap_", var, "_", ref, ".jpeg")),
        width = 10, height = 9, units = "in", res = 600, quality = 95)
    layout(matrix(1:2, 2, 1), heights = c(1, 0.06))
    par(mar = c(1, 2, 3, 2))
    plot(r_ov, col = col_ov, breaks = ov_brks,
         main = "",
         axes = FALSE, legend = FALSE, colNA = "white")
    draw_borders(ref)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, cex = 1.0, bty = "n", horiz = TRUE)
    dev.off()
  }

  ov_df <- do.call(rbind, overlap_rows)
  write.csv(ov_df, file.path(config$plot_dir,
                             paste0(id, "_overlap_summary_", ref, ".csv")), row.names = FALSE)

  message(sprintf("\n    OVERLAP SUMMARY [%s, %s]:", ref, thresholds[[main_threshold]]$label))
  message(sprintf("    %-30s %8s %12s %8s", "Variable", "Cells", "PLHIV", "%SSA"))
  for (i in seq_len(nrow(ov_df))) {
    message(sprintf("    %-30s %8s %12s %7.1f%%",
                    ov_df$label[i],
                    format(ov_df$n_both_hotspot[i], big.mark = ","),
                    format(ov_df$plhiv_both_hotspot[i], big.mark = ","),
                    ifelse(is.na(ov_df$pct_plhiv_ssa[i]), 0, ov_df$pct_plhiv_ssa[i])))
  }


  # ============================================================
  # 2. PANEL MAPS + STANDALONE BH MAPS
  # ============================================================

  # --- Panel A: 3 main EWEs at 95%   1 row, labelled (a)/(b)/(c) ---
  main_ov_rasters <- list()
  for (var in config$main_vars) {
    ewe_nm <- cat_layer_name(var)
    if (!ewe_nm %in% names(r_cat)) next
    ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
    ov <- make_overlap(ewe_v, hiv_v, thresholds[[main_threshold]]$cat_min)
    r_ov <- rast(hiv_cat); values(r_ov) <- ov
    main_ov_rasters[[var]] <- r_ov
  }

  if (length(main_ov_rasters) >= 2) {
    nv <- length(main_ov_rasters)
    jpeg(file.path(config$plot_dir, "main",
                   paste0(id, "_panel_main3_95pct_", ref, ".jpeg")),
         width = 5 * nv, height = 6, units = "in", res = 600, quality = 95)
    layout(rbind(seq_len(nv), rep(nv + 1, nv)), heights = c(1, 0.12))
    par(mar = c(1, 1, 3, 1), oma = c(0, 0, 3, 0))
    abc <- letters[seq_len(nv)]
    for (i in seq_along(main_ov_rasters)) {
      var <- names(main_ov_rasters)[i]
      plot(main_ov_rasters[[i]], col = col_ov, breaks = ov_brks,
           main = paste0(abc[i], ") ", get_label(names(main_ov_rasters)[i])),
           axes = FALSE, legend = FALSE, colNA = "white", cex.main = 1.1)
      draw_borders(ref)
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, ncol = 2, cex = 1.0, bty = "n",
           x.intersp = 0.5, text.width = strwidth("Not an overlapping hotspot") * 1.2)
    dev.off()
    message("    3-main panel (a/b/c) at 95% saved")
  }

  # --- Panel B: Multi-threshold   ROWS = EWE, COLS = threshold (horizontal) ---
  if (length(config$main_vars) >= 2) {
    nv <- length(config$main_vars)
    nthr <- 3
    jpeg(file.path(config$plot_dir, "main",
                   paste0(id, "_panel_main_multithreshold_", ref, ".jpeg")),
         width = 5 * nthr, height = 4 * nv + 3, units = "in", res = 600, quality = 95)
    layout(rbind(matrix(seq_len(nv * nthr), nv, nthr, byrow = TRUE), rep(nv * nthr + 1, nthr)), heights = c(rep(1, nv), 0.12))
    par(mar = c(1, 1, 3, 1), oma = c(0, 0, 4, 0))

    abc <- letters
    vi <- 1
    for (var in config$main_vars) {
      ewe_nm <- cat_layer_name(var)
      for (thr_name in c("90", "95", "99")) {
        thr <- thresholds[[thr_name]]
        if (!ewe_nm %in% names(r_cat)) { plot.new(); vi <- vi + 1; next }
        ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
        ov <- make_overlap(ewe_v, hiv_v, thr$cat_min)
        r_ov <- rast(hiv_cat); values(r_ov) <- ov
        plot(r_ov, col = col_ov, breaks = ov_brks,
             main = paste0(abc[vi], ") ", get_label(var), " (Gi* ", thr$label, ")"),
             axes = FALSE, legend = FALSE, colNA = "white", cex.main = 0.9)
        draw_borders(ref)
        vi <- vi + 1
      }
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, ncol = 2, cex = 1.0, bty = "n",
           x.intersp = 0.5, text.width = strwidth("Not an overlapping hotspot") * 1.2)
    dev.off()
    message("    Multi-threshold panel (rows=EWE, cols=threshold) saved")
  }

  # --- Panel C: Sensitivity DROUGHT  x  3 thresholds ---
  drought_sens <- c("SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0")
  drought_present <- present_cat_vars(drought_sens)

  if (length(drought_present) >= 2) {
    jpeg(file.path(config$plot_dir, "robustness",
                   paste0(id, "_panel_drought_sens_multithr_", ref, ".jpeg")),
         width = 15, height = 4 * length(drought_present) + 3, units = "in", res = 600, quality = 95)
    nd <- length(drought_present)
    layout(rbind(matrix(seq_len(nd * 3), nd, 3, byrow = TRUE), rep(nd * 3 + 1, 3)), heights = c(rep(1, nd), 0.12))
    par(mar = c(1, 1, 3, 1), oma = c(0, 0, 4, 0))
    abc <- letters; vi <- 1
    for (var in drought_present) {
      ewe_nm <- cat_layer_name(var)
      for (thr_name in c("90", "95", "99")) {
        thr <- thresholds[[thr_name]]
        if (!ewe_nm %in% names(r_cat)) { plot.new(); vi <- vi + 1; next }
        ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
        ov <- make_overlap(ewe_v, hiv_v, thr$cat_min)
        r_ov <- rast(hiv_cat); values(r_ov) <- ov
        plot(r_ov, col = col_ov, breaks = ov_brks,
             main = paste0(abc[vi], ") ", get_label(var), " (Gi* ", thr$label, ")"),
             axes = FALSE, legend = FALSE, colNA = "white", cex.main = 0.9)
        draw_borders(ref)
        vi <- vi + 1
      }
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, ncol = 2, cex = 1.0, bty = "n",
           x.intersp = 0.5, text.width = strwidth("Not an overlapping hotspot") * 1.2)
    dev.off()
    message("    Drought sensitivity multi-threshold panel saved")
  }

  # --- Panel D: Other sensitivity (SPI-1.5, TX95p, R99p)  x  3 thresholds ---
  other_sens <- c("TX95p_HWMF", "R99pTOT", "SPI3_leq1p5")
  other_present <- present_cat_vars(other_sens)

  if (length(other_present) >= 2) {
    jpeg(file.path(config$plot_dir, "robustness",
                   paste0(id, "_panel_other_sens_multithr_", ref, ".jpeg")),
         width = 15, height = 4 * length(other_present) + 3, units = "in", res = 600, quality = 95)
    no <- length(other_present)
    layout(rbind(matrix(seq_len(no * 3), no, 3, byrow = TRUE), rep(no * 3 + 1, 3)), heights = c(rep(1, no), 0.12))
    par(mar = c(1, 1, 3, 1), oma = c(0, 0, 4, 0))
    abc <- letters; vi <- 1
    for (var in other_present) {
      ewe_nm <- cat_layer_name(var)
      for (thr_name in c("90", "95", "99")) {
        thr <- thresholds[[thr_name]]
        if (!ewe_nm %in% names(r_cat)) { plot.new(); vi <- vi + 1; next }
        ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
        ov <- make_overlap(ewe_v, hiv_v, thr$cat_min)
        r_ov <- rast(hiv_cat); values(r_ov) <- ov
        plot(r_ov, col = col_ov, breaks = ov_brks,
             main = paste0(abc[vi], ") ", get_label(var), " (Gi* ", thr$label, ")"),
             axes = FALSE, legend = FALSE, colNA = "white", cex.main = 0.9)
        draw_borders(ref)
        vi <- vi + 1
      }
    }
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, ncol = 2, cex = 1.0, bty = "n",
           x.intersp = 0.5, text.width = strwidth("Not an overlapping hotspot") * 1.2)
    dev.off()
    message("    Other sensitivity multi-threshold panel saved")
  }

  # --- Panel E: All 6 sensitivity EWEs at 95% (single threshold) ---
  all_sens <- c("TX95p_HWMF", "R99pTOT", "SPI3_leq1p5",
                "SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0")
  sens_ov_rasters <- list()
  for (v in all_sens) {
    ewe_nm <- cat_layer_name(v)
    if (!ewe_nm %in% names(r_cat)) next
    ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
    ov <- make_overlap(ewe_v, hiv_v, thresholds[[main_threshold]]$cat_min)
    if (sum(ov == 1, na.rm = TRUE) > 0) {
      r_ov <- rast(hiv_cat); values(r_ov) <- ov
      sens_ov_rasters[[v]] <- r_ov
    }
  }
  if (length(sens_ov_rasters) >= 2) {
    ns <- length(sens_ov_rasters)
    nr_s <- ceiling(ns / 3); nc_s <- min(3, ns)
    jpeg(file.path(config$plot_dir, "robustness",
                   paste0(id, "_panel_all6_sens_95pct_", ref, ".jpeg")),
         width = 5 * nc_s, height = 4 * nr_s + 3, units = "in", res = 600, quality = 95)
    layout(rbind(matrix(seq_len(nr_s * nc_s), nr_s, nc_s, byrow = TRUE),
                 rep(nr_s * nc_s + 1, nc_s)), heights = c(rep(1, nr_s), 0.08))
    par(mar = c(1, 1, 2.5, 1), oma = c(0, 0, 0, 0))
    abc_s <- letters[seq_along(sens_ov_rasters)]
    for (si in seq_along(sens_ov_rasters)) {
      plot(sens_ov_rasters[[si]], col = col_ov, breaks = ov_brks,
           main = paste0(abc_s[si], ") ", get_label(names(sens_ov_rasters)[si])),
           axes = FALSE, legend = FALSE, colNA = "white", cex.main = 1.0)
      draw_borders(ref)
    }
    remainder_s <- (nc_s * nr_s) - ns
    if (remainder_s > 0) for (ri in seq_len(remainder_s)) plot.new()
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, ncol = 2, cex = 1.0, bty = "n",
           x.intersp = 0.5)
    dev.off()
    message("    All 6-sensitivity panel [", ref, "] saved")
  }

  # --- Standalone BH maps for each variable ---
  for (var in all_vars) {
    ewe_nm <- cat_layer_name(var)
    if (!ewe_nm %in% names(r_cat)) next
    ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)
    ov <- make_overlap(ewe_v, hiv_v, thresholds[[main_threshold]]$cat_min)
    r_ov <- rast(hiv_cat); values(r_ov) <- ov
    is_main <- var %in% config$main_vars
    subdir <- if (is_main) "main" else "robustness"

    jpeg(file.path(config$plot_dir, subdir,
                   paste0(id, "_standalone_", var, "_", ref, ".jpeg")),
         width = 10, height = 9, units = "in", res = 600, quality = 95)
    layout(matrix(1:2, 2, 1), heights = c(1, 0.06))
    par(mar = c(1, 2, 3, 2))
    plot(r_ov, col = col_ov, breaks = ov_brks,
         main = "",
         axes = FALSE, legend = FALSE, colNA = "white")
    draw_borders(ref)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center", legend = ov_leg, fill = col_ov, cex = 1.0, bty = "n", horiz = TRUE)
    dev.off()
  }
  message("    Standalone BH maps saved")


  # ============================================================
  # 3. PLHIV at every threshold, continent-wide and per country, with intervals
  # ============================================================

  # Country code -> name, read from the country layer of the draws stack. Stops
  # when the table is absent rather than labelling countries by number.
  if (!file.exists(config$plhiv_draws_path))
    stop("The draws stack is missing, so countries cannot be named: ",
         config$plhiv_draws_path)
  tmp_r <- rast(config$plhiv_draws_path)
  cname_levels <- get_country_levels(find_country_layer(tmp_r))
  rm(tmp_r)
  message(sprintf("    Country lookup: %d countries", nrow(cname_levels)))

  get_cname <- function(cid) country_name(cid, cname_levels)

  plhiv_thr_rows <- list()

  for (var in all_vars) {
    ewe_nm <- cat_layer_name(var)
    if (!ewe_nm %in% names(r_cat)) next
    ewe_v <- values(r_cat[[ewe_nm]], mat = FALSE)

    for (thr_name in names(thresholds)) {
      thr <- thresholds[[thr_name]]
      ov <- make_overlap(ewe_v, hiv_v, thr$cat_min)
      both_mask <- !is.na(ov) & ov == 1

      # SSA total with CI
      ci_ssa <- count_plhiv_ci(both_mask)
      pci_ssa <- pct_plhiv_ci(both_mask, !is.na(hiv_v))
      plhiv_thr_rows[[length(plhiv_thr_rows) + 1]] <- data.frame(
        mode = ref, variable = var, label = get_label(var),
        analysis = if (var %in% config$main_vars) "main" else "robustness",
        threshold = thr$label, region = "SSA",
        n_both_hotspot = sum(both_mask, na.rm = TRUE),
        plhiv_mean = ci_ssa$plhiv_mean,
        plhiv_lower = ci_ssa$plhiv_lower,
        plhiv_upper = ci_ssa$plhiv_upper,
        pct_mean = pci_ssa$pct_mean,
        pct_lower = pci_ssa$pct_lower,
        pct_upper = pci_ssa$pct_upper,
        stringsAsFactors = FALSE)

      # Per country with CI
      if (!is.null(country_v)) {
        cids <- sort(unique(country_v[!is.na(country_v) & country_v > 0]))
        for (cid in cids) {
          c_mask <- !is.na(country_v) & country_v == cid
          c_both <- c_mask & both_mask
          ci_c <- count_plhiv_ci(c_both)
          ci_c_total <- count_plhiv_ci(c_mask)
          pct_nat <- if (ci_c_total$plhiv_mean > 0) round(100 * ci_c$plhiv_mean / ci_c_total$plhiv_mean, 2) else NA

          pci_c <- pct_plhiv_ci(c_both, c_mask)
          plhiv_thr_rows[[length(plhiv_thr_rows) + 1]] <- data.frame(
            mode = ref, variable = var, label = get_label(var),
            analysis = if (var %in% config$main_vars) "main" else "robustness",
            threshold = thr$label, region = get_cname(cid),
            n_both_hotspot = sum(c_both, na.rm = TRUE),
            plhiv_mean = ci_c$plhiv_mean,
            plhiv_lower = ci_c$plhiv_lower,
            plhiv_upper = ci_c$plhiv_upper,
            pct_mean = pci_c$pct_mean,
            pct_lower = pci_c$pct_lower,
            pct_upper = pci_c$pct_upper,
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  plhiv_thr_df <- do.call(rbind, plhiv_thr_rows)
  write.csv(plhiv_thr_df, file.path(config$plot_dir,
            paste0(id, "_overlap_plhiv_all_thresholds_", ref, ".csv")), row.names = FALSE)

  ssa_sub <- plhiv_thr_df[plhiv_thr_df$region == "SSA", ]
  message(sprintf("\n    PLHIV IN HOTSPOT OVERLAPS   SSA [%s] (mean [95%% UI]):", ref))
  for (i in seq_len(nrow(ssa_sub))) {
    lo <- ssa_sub$plhiv_lower[i]; hi <- ssa_sub$plhiv_upper[i]
    ci_txt <- if (!is.na(lo)) sprintf(" [%s\u2013%s]", format(lo, big.mark=","), format(hi, big.mark=",")) else ""
    message(sprintf("    %-30s %6s %8s %12s%s",
                    ssa_sub$label[i], ssa_sub$threshold[i],
                    format(ssa_sub$n_both_hotspot[i], big.mark = ","),
                    format(ssa_sub$plhiv_mean[i], big.mark = ","), ci_txt))
  }

  rm(r_cat, r_res); gc(FALSE)
}

# ============================================================
# RUN
# ============================================================

if (length(cat_ssa) && length(res_ssa))
  process_mode(cat_ssa[1], res_ssa[1], "SSA")

if (length(cat_cty) && length(res_cty))
  process_mode(cat_cty[1], res_cty[1], "COUNTRY")


message("Done. Pipeline complete.")
