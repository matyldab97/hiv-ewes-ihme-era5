# ============================================================
# Indicator screening
# ------------------------------------------------------------
# Maps, distributions and severity classifications for every candidate
# indicator, used to select the indicators reported in the main analysis.
#
# Requires: steps 02 to 12
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))
source(file.path("R", "01_setup", "03_lib_shared.R"))

# One page per indicator, held in EWE_verification_figures as a PDF and a JPEG:
# map, histogram, interquartile classification for the whole of SSA, and the
# same classification applied within each country.
#
# Class boundaries fall into the higher class. Severity sums are
# negated so that higher always means more severe. A cell is assigned to
# the country covering the largest share of its area.
#
# The workbook holds continent-wide cell counts on the first sheet and
# per-country counts on the second.

library(terra)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(openxlsx)
library(sf)
library(exactextractr)

apply_mask <- function(r, mask_r) {
  mask_res <- resample(mask_r, r, method = "near")
  r_masked <- mask(r, mask_res)
  inside <- !is.na(mask_res)
  r_masked[inside & is.na(r_masked)] <- 0
  r_masked
}

raster_to_png_grob <- function(r, cols, main_title, breaks = NULL,
                               cex_main = 0.8, borders = NULL, cat_labels = NULL) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = 1320, height = 1000)
  if (!is.null(cat_labels)) {
    rr <- r
    levels(rr) <- data.frame(id = seq_along(cat_labels) - 1L, class = cat_labels)
    par(mar = c(1, 1, 2, 3))
    plot(rr, col = cols, main = main_title, axes = FALSE, box = FALSE,
         cex.main = cex_main, mar = c(1, 1, 2, 3),
         plg = list(cex = 2.8, x = "bottom", ncol = 5, bty = "n", x.intersp = 0.4))
    if (!is.null(borders)) lines(borders, col = "black", lwd = 0.6)
  } else if (!is.null(breaks)) {
    par(mar = c(1, 1, 2, 4))
    plot(r, col = cols, breaks = breaks, main = main_title,
         axes = FALSE, box = FALSE, cex.main = cex_main)
    if (!is.null(borders)) lines(borders, col = "black", lwd = 0.6)
  } else {
    par(mar = c(1, 1, 2, 4))
    plot(r, col = cols, main = main_title,
         axes = FALSE, box = FALSE, cex.main = cex_main)
    if (!is.null(borders)) lines(borders, col = "black", lwd = 0.6)
  }
  dev.off()
  img <- png::readPNG(tmp)
  rasterGrob(img, interpolate = TRUE)
}

make_labels <- function(br, k_eff, prefix = "C", discrete = FALSE) {
  labs <- character(k_eff)
  for (i in seq_len(k_eff)) {
    if (discrete) {
      lo <- if (i == 1) ceiling(br[1]) else as.integer(ceiling(br[i]))
      if (i < k_eff) {
        hi_raw <- br[i + 1]
        hi <- if (hi_raw == floor(hi_raw)) as.integer(hi_raw) - 1L else floor(hi_raw)
      } else { hi <- floor(br[i + 1]) }
      if (is.na(lo) || is.na(hi)) { labs[i] <- paste0(prefix, i); next }
      if (lo > hi) hi <- lo
      if (lo == hi) labs[i] <- paste0(prefix, i, ": ", lo)
      else labs[i] <- paste0(prefix, i, ": ", lo, "-", hi)
    } else {
      labs[i] <- paste0(prefix, i, ": ", fmt3(br[i]), "-", fmt3(br[i + 1]))
    }
  }
  labs
}

# --- Country raster: majority rule ---
count_classes <- function(cls, max_k = 4) {
  counts <- rep(NA_integer_, max_k)
  if (length(cls) == 0) return(counts)
  for (i in seq_len(max_k)) counts[i] <- sum(cls == i, na.rm = TRUE)
  counts
}

# --- Panel functions ---
panel_raw_map <- function(r_masked, nm, cex_main = 1.0) {
  raster_to_png_grob(r_masked, c("grey90", rev(heat.colors(100))), "", cex_main = cex_main)
}

panel_histogram <- function(r_masked, nm) {
  vals <- as.vector(values(r_masked, na.rm = TRUE))
  vals_pos <- vals[vals != 0]
  if (length(vals_pos) == 0) return(textGrob(paste0(nm, "\nNo exposed cells")))
  q1 <- quantile(vals_pos, 0.25); med <- quantile(vals_pos, 0.50); q3 <- quantile(vals_pos, 0.75)
  df <- data.frame(value = vals_pos)
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = q1,  linetype = "dashed", colour = "orange", linewidth = 0.8) +
    geom_vline(xintercept = med, linetype = "solid",  colour = "red",    linewidth = 1.0) +
    geom_vline(xintercept = q3,  linetype = "dashed", colour = "orange", linewidth = 0.8) +
    labs(title = NULL, x = "Value", y = "Frequency") + theme_minimal() +
    theme(plot.margin = margin(t = 28, r = 45, b = 6, l = 28))
  ggplotGrob(p)
}

panel_iqr_classes <- function(r_masked, nm) {
  v_all <- as.vector(values(r_masked, na.rm = FALSE))
  v_pos <- v_all[!is.na(v_all) & v_all != 0]
  rc <- r_masked
  if (length(v_pos) == 0) { rc[] <- NA; return(raster_to_png_grob(rc, "grey", paste0(nm, " - IQR/Median (no data)"))) }
  br <- compute_iqr_breaks(v_pos); k_eff <- length(br) - 1
  class_ids <- rep(NA_integer_, length(v_all)); class_ids[!is.na(v_all) & v_all == 0] <- 0L
  idx <- !is.na(v_all) & v_all != 0
  if (k_eff >= 1 && any(idx)) class_ids[idx] <- classify_vals(v_all[idx], br)
  else if (any(idx)) { class_ids[idx] <- 1L; k_eff <- 1 }
  values(rc) <- class_ids
  cols <- c("grey85", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d")[1:(k_eff + 1)]
  sev_lab <- c("Unexposed", "Low", "Moderate", "High", "Very High")[1:(k_eff + 1)]
  raster_to_png_grob(rc, cols, "", cat_labels = sev_lab)
}


panel_country_classes <- function(r_masked, nm, country_assign, borders_v, method = "iqr") {
  v_all <- as.vector(values(r_masked, na.rm = FALSE)); n <- length(v_all)
  class_ids <- rep(NA_integer_, n); class_ids[!is.na(v_all) & v_all == 0] <- 0L
  countries <- sort(unique(na.omit(country_assign)))
  for (cn in countries) {
    cidx <- which(country_assign == cn & !is.na(v_all) & v_all != 0)
    if (length(cidx) == 0) next
    vp <- v_all[cidx]
    br <- if (method == "iqr") compute_iqr_breaks(vp) else compute_range_breaks(vp)
    if (is.null(br) || length(br) < 2) next
    class_ids[cidx] <- classify_vals(vp, br)
  }
  rc <- r_masked; values(rc) <- class_ids
  if (method == "iqr") {
    title_str <- ""
    cols <- c("grey85", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d")
  } else {
    title_str <- ""
    cols <- c("grey85", "#eff3ff", "#6baed6", "#2171b5", "#08306b")
  }
  present <- sort(unique(na.omit(class_ids))); present <- present[present > 0]
  max_cls <- if (length(present) > 0) max(present) else 4
  n_cls <- min(max_cls, 4)
  sev_lab <- c("Unexposed", "Low", "Moderate", "High", "Very High")[1:(n_cls + 1)]
  raster_to_png_grob(rc, cols[1:(n_cls + 1)], title_str, borders = borders_v, cat_labels = sev_lab)
}

# --- Stats helpers ---
compute_global_stats <- function(vp) {
  iq_br <- compute_iqr_breaks(vp); rg_br <- compute_range_breaks(vp); discrete <- is_discrete(vp)
  iq_counts <- rep(NA_integer_, 4); iq_str <- ""
  if (!is.null(iq_br) && length(iq_br) >= 2) {
    k <- length(iq_br) - 1; cls <- classify_vals(vp, iq_br); iq_counts <- count_classes(cls)
    iq_str <- paste(make_labels(iq_br, k, "C", discrete), collapse = " | ")
  }
  rg_counts <- rep(NA_integer_, 4); rg_str <- ""
  if (!is.null(rg_br) && length(rg_br) >= 2) {
    k <- length(rg_br) - 1; cls <- classify_vals(vp, rg_br); rg_counts <- count_classes(cls)
    rg_str <- paste(make_labels(rg_br, k, "R", discrete), collapse = " | ")
  }
  list(iq_counts = iq_counts, iq_str = iq_str, rg_counts = rg_counts, rg_str = rg_str)
}

compute_country_stats <- function(v_all, country_assign, nm) {
  countries <- sort(unique(na.omit(country_assign))); rows <- list()
  for (cn in countries) {
    cidx <- which(country_assign == cn & !is.na(v_all))
    if (length(cidx) == 0) next
    cv <- v_all[cidx]; n_total <- length(cv); n_exposed <- sum(cv != 0); n_unexposed <- sum(cv == 0)
    vp <- cv[cv != 0]; gs <- compute_global_stats(vp)
    row <- data.frame(
      Indicator = nm, Country = cn, Total_Cells = n_total,
      Exposed = n_exposed, Unexposed = n_unexposed,
      Pct_Exposed = round(100 * n_exposed / max(1, n_total), 1),
      Min = if (length(vp) > 0) round(min(vp), 3) else NA,
      Q1 = if (length(vp) > 0) round(quantile(vp, 0.25), 3) else NA,
      Median = if (length(vp) > 0) round(quantile(vp, 0.50), 3) else NA,
      Q3 = if (length(vp) > 0) round(quantile(vp, 0.75), 3) else NA,
      Max = if (length(vp) > 0) round(max(vp), 3) else NA,
      IQR_C1_n = gs$iq_counts[1], IQR_C2_n = gs$iq_counts[2],
      IQR_C3_n = gs$iq_counts[3], IQR_C4_n = gs$iq_counts[4], IQR_Breaks = gs$iq_str,
      Range_R1_n = gs$rg_counts[1], Range_R2_n = gs$rg_counts[2],
      Range_R3_n = gs$rg_counts[3], Range_R4_n = gs$rg_counts[4], Range_Breaks = gs$rg_str,
      stringsAsFactors = FALSE)
    rows <- append(rows, list(row))
  }
  if (length(rows) > 0) do.call(rbind, rows) else NULL
}

# --- Metric pages: page 1 global + page 2 country ---
metric_pages <- function(r_masked, nm, country_assign, borders_v) {
  ## Four panels: a) raw map, b) histogram, c) SSA-wide d) within-country
  g_raw  <- panel_raw_map(r_masked, nm)
  g_hist <- panel_histogram(r_masked, nm)
  g_ssa  <- panel_iqr_classes(r_masked, nm)
  v_all <- as.vector(values(r_masked, na.rm = FALSE))
  g_last <- panel_country_classes(r_masked, nm, country_assign, borders_v, "iqr")
  grobs <- list(g_raw, g_hist, g_ssa, g_last); nc <- 2; nr <- 2
  ## helper: draw the 2x2 grid + a)/b)/c)/d) corner labels, no titles
  draw_labelled <- function() {
    do.call(grid.arrange, c(grobs, list(ncol = nc, nrow = nr)))
    labs <- c("a)", "b)", "c)", "d)")
    xs <- c(0.010, 0.505, 0.010, 0.505); ys <- c(0.992, 0.992, 0.498, 0.498)
    for (k in seq_along(labs))
      grid.text(labs[k], x = xs[k], y = ys[k], just = c("left","top"),
                gp = gpar(fontsize = 13, fontface = "bold"))
  }
  grid.newpage(); draw_labelled()
  ## individual single-page figure (for the supplement)
  if (exists("EWE_FIG_DIR")) {
    safe <- gsub("[^A-Za-z0-9]+", "_", nm)
    jpeg(file.path(EWE_FIG_DIR, paste0(safe, ".jpeg")),
         width = 2700, height = 2050, res = 300, quality = 95)
    draw_labelled(); dev.off()
    ## individual PDF (same square layout) - upload THIS, no pixel limit
    pdf(file.path(EWE_FIG_DIR, paste0(safe, ".pdf")), width = 10.8, height = 8.2)
    draw_labelled(); dev.off()
  }
  vals <- as.vector(values(r_masked, na.rm = TRUE)); vp <- vals[vals != 0]
  n_exposed <- length(vp); n_unexposed <- sum(vals == 0); n_total <- n_exposed + n_unexposed
  gs <- compute_global_stats(vp)
  global_row <- data.frame(
    Indicator = nm, Total_Cells = n_total, Exposed = n_exposed, Unexposed = n_unexposed,
    Pct_Exposed = round(100 * n_exposed / max(1, n_total), 1),
    Min = if (length(vp) > 0) round(min(vp), 3) else NA,
    Q1 = if (length(vp) > 0) round(quantile(vp, 0.25), 3) else NA,
    Median = if (length(vp) > 0) round(quantile(vp, 0.50), 3) else NA,
    Q3 = if (length(vp) > 0) round(quantile(vp, 0.75), 3) else NA,
    Max = if (length(vp) > 0) round(max(vp), 3) else NA,
    IQR_C1_n = gs$iq_counts[1], IQR_C2_n = gs$iq_counts[2],
    IQR_C3_n = gs$iq_counts[3], IQR_C4_n = gs$iq_counts[4], IQR_Breaks = gs$iq_str,
    Range_R1_n = gs$rg_counts[1], Range_R2_n = gs$rg_counts[2],
    Range_R3_n = gs$rg_counts[3], Range_R4_n = gs$rg_counts[4], Range_Breaks = gs$rg_str,
    stringsAsFactors = FALSE)
  country_detail <- compute_country_stats(v_all, country_assign, nm)
  list(global = global_row, country = country_detail)
}

# --- Comparison pages ---
# --- Main engine ---
run_full_analysis <- function(metric_files, mask_file, xlsx_out,
                              precomputed = NULL, country_assign, borders_v) {
  .mstk <- rast(mask_file)
  .pi <- grep("^PLHIV_mean$", names(.mstk))
  if (length(.pi) == 0) .pi <- grep("^PLHIV", names(.mstk))
  if (length(.pi) == 0) stop("The interval raster has no PLHIV layer, so the analysis domain cannot be derived; rerun step 12.")
  mask_r <- .mstk[[.pi[1]]]; stored_rasters <- list()
  global_stats <- list(); country_stats <- list()
  ## folder for the individual single-page figures, one PDF and one JPEG per
  ## indicator, written inside metric_pages
  EWE_FIG_DIR <<- file.path(DIR_SCREENING, "EWE_verification_figures")
  dir.create(EWE_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  pdf(NULL, width = 10.8, height = 8.2)
  for (nm in names(metric_files)) {
    cat("Processing:", nm, "...")
    r <- rast(metric_files[[nm]])
    if (grepl("Severity", nm, ignore.case = TRUE)) { r <- r * -1; cat(" [flipped]") }
    r_masked <- apply_mask(r, mask_r); stored_rasters[[nm]] <- r_masked
    vals <- as.vector(values(r_masked, na.rm = TRUE))
    cat(" cells:", length(vals), " exposed:", sum(vals != 0), "\n")
    res <- metric_pages(r_masked, nm, country_assign, borders_v)
    global_stats <- append(global_stats, list(res$global))
    if (!is.null(res$country)) country_stats <- append(country_stats, list(res$country))
  }
  if (!is.null(precomputed)) {
    for (nm in names(precomputed)) {
      cat("Processing (precomputed):", nm, "...")
      r <- precomputed[[nm]]
      if (grepl("Severity", nm, ignore.case = TRUE)) { r <- r * -1; cat(" [flipped]") }
      r_masked <- apply_mask(r, mask_r); stored_rasters[[nm]] <- r_masked
      vals <- as.vector(values(r_masked, na.rm = TRUE))
      cat(" cells:", length(vals), " exposed:", sum(vals != 0), "\n")
      res <- metric_pages(r_masked, nm, country_assign, borders_v)
      global_stats <- append(global_stats, list(res$global))
      if (!is.null(res$country)) country_stats <- append(country_stats, list(res$country))
    }
  }
  dev.off(); cat("Per-indicator figures:", EWE_FIG_DIR, "\n")
  wb <- createWorkbook()
  global_df <- do.call(rbind, global_stats)
  addWorksheet(wb, "Global"); writeData(wb, "Global", global_df)
  setColWidths(wb, "Global", cols = 1:ncol(global_df), widths = "auto")
  freezePane(wb, "Global", firstRow = TRUE)
  if (length(country_stats) > 0) {
    country_df <- do.call(rbind, country_stats)
    addWorksheet(wb, "Country_Detail"); writeData(wb, "Country_Detail", country_df)
    setColWidths(wb, "Country_Detail", cols = 1:ncol(country_df), widths = "auto")
    freezePane(wb, "Country_Detail", firstRow = TRUE)
  }
  saveWorkbook(wb, xlsx_out, overwrite = TRUE); cat("Excel saved:", xlsx_out, "\n")
  invisible(global_df)
}

# ============================================================
# SETUP: COUNTRY LAYER (computed once, reused for all 4 PDFs)
# ============================================================
mask_file <- file.path(DIR_MC, sprintf("ALL_15_59_%d_UI_0p25.tif", TARGET_YEAR))  # analysis domain, taken from the modelled PLHIV footprint
gadm_path <- FILE_GADM

cat("Building country assignment layer...\n")
ewe_ref <- rast(file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWF/TX90p_ANN_HWF_%d.tif", TARGET_YEAR)))
tpl <- ewe_ref; tpl[] <- 1:ncell(tpl)
country_info <- compute_country_raster(tpl, gadm_path)
country_assign <- country_info$assignment
borders_v <- vect(country_info$countries_sf)
cat("Countries found:", length(unique(na.omit(country_assign))), "\n")
cat("Cells assigned:", sum(!is.na(country_assign)), "/", length(country_assign), "\n")

# ============================================================
# Drought, 32 indicators
# ============================================================
drought_metrics <- list(
  "SPI3 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SPI1 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SPI/SPI1_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SPEI3 Month Count (leq-1.0)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SPEI1 Month Count (leq-1.0)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI1_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SRI3 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SRI1 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SRI/SRI1_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SMA3 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SMA1 Month Count (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SMA/SMA1_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
  "SPI3 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SPI1 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SPI/SPI1_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SPEI3 Severity Sum (leq-1.0)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SPEI1 Severity Sum (leq-1.0)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI1_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SRI3 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SRI1 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SRI/SRI1_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SMA3 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SMA1 Severity Sum (leq-1.0)"  = file.path(DIR_INDICES, sprintf("SMA/SMA1_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
  "SPI3 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SPI1 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SPI/SPI1_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SPEI3 Month Count (leq-1.5)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SPEI1 Month Count (leq-1.5)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI1_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SRI3 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SRI1 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SRI/SRI1_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SMA3 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SMA1 Month Count (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SMA/SMA1_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
  "SPI3 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SPI1 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SPI/SPI1_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SPEI3 Severity Sum (leq-1.5)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SPEI1 Severity Sum (leq-1.5)" = file.path(DIR_INDICES, sprintf("SPEI/SPEI1_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SRI3 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SRI1 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SRI/SRI1_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SMA3 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
  "SMA1 Severity Sum (leq-1.5)"  = file.path(DIR_INDICES, sprintf("SMA/SMA1_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR))
)

run_full_analysis(drought_metrics, mask_file,
                  file.path(DIR_SCREENING, "Drought_Indicators.xlsx"),
                  country_assign = country_assign, borders_v = borders_v)

# ============================================================
# Heatwaves, 14 indicators
# ============================================================
heatwave_metrics <- list(
  # TX90p metrics
  "TX90p HWF (heatwave days)"       = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWF/TX90p_ANN_HWF_%d.tif", TARGET_YEAR)),
  "TX90p HWN (events)"              = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWN/TX90p_ANN_HWN_%d.tif", TARGET_YEAR)),
  "TX90p HWM (mean magnitude C)"    = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWM/TX90p_ANN_HWM_%d.tif", TARGET_YEAR)),
  "TX90p HWMF (magnitude x days)"   = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWMF/TX90p_ANN_HWMF_%d.tif", TARGET_YEAR)),
  "TX90p HWNM (events x magnitude)" = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWNM/TX90p_ANN_HWNM_%d.tif", TARGET_YEAR)),
  "TX90p HWES (excess heat C-days)"  = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWES/TX90p_ANN_HWES_%d.tif", TARGET_YEAR)),
  "TX90p HWPD (pctile deviation)"    = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWPD/TX90p_ANN_HWPD_%d.tif", TARGET_YEAR)),
  # TX95p metrics
  "TX95p HWF (heatwave days)"       = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWF/TX95p_ANN_HWF_%d.tif", TARGET_YEAR)),
  "TX95p HWN (events)"              = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWN/TX95p_ANN_HWN_%d.tif", TARGET_YEAR)),
  "TX95p HWM (mean magnitude C)"    = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWM/TX95p_ANN_HWM_%d.tif", TARGET_YEAR)),
  "TX95p HWMF (magnitude x days)"   = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWMF/TX95p_ANN_HWMF_%d.tif", TARGET_YEAR)),
  "TX95p HWNM (events x magnitude)" = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWNM/TX95p_ANN_HWNM_%d.tif", TARGET_YEAR)),
  "TX95p HWES (excess heat C-days)"  = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWES/TX95p_ANN_HWES_%d.tif", TARGET_YEAR)),
  "TX95p HWPD (pctile deviation)"    = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWPD/TX95p_ANN_HWPD_%d.tif", TARGET_YEAR))
)

run_full_analysis(heatwave_metrics, mask_file,
                  file.path(DIR_SCREENING, "Heatwave_Indicators.xlsx"),
                  country_assign = country_assign, borders_v = borders_v)

# ============================================================
# Extreme rainfall, 8 indicators
# ============================================================
rainfall_file_metrics <- list(
  # R95p metrics
  "R95pDAY (extreme days)"            = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_DAY/R95pDAY_%d.tif", TARGET_YEAR)),
  "R95pTOT (total extreme rain)"      = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_TOT/R95pTOT_%d.tif", TARGET_YEAR)),
  "R95p Excess Precip (mm)"           = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_EXCESS/R95pEXCESS_%d.tif", TARGET_YEAR)),
  "R95p Pctile Deviation (pctile-days)" = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_PCTLDEV/R95pPCTD_%d.tif", TARGET_YEAR)),
  # R99p metrics
  "R99pDAY (extreme days)"            = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_DAY/R99pDAY_%d.tif", TARGET_YEAR)),
  "R99pTOT (total extreme rain)"      = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_TOT/R99pTOT_%d.tif", TARGET_YEAR)),
  "R99p Excess Precip (mm)"           = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_EXCESS/R99pEXCESS_%d.tif", TARGET_YEAR)),
  "R99p Pctile Deviation (pctile-days)" = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_PCTLDEV/R99pPCTD_%d.tif", TARGET_YEAR))
)

run_full_analysis(rainfall_file_metrics, mask_file,
                  file.path(DIR_SCREENING, "Rainfall_Indicators.xlsx"),
                  country_assign = country_assign, borders_v = borders_v)

message("Done. Next: step 14.")
