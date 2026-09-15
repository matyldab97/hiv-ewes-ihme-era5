# ============================================================
# PLHIV exposure classification
# ------------------------------------------------------------
# Part 1  single hazard, 12 indicators, four classification scenarios
# Part 2  figures and tables
#
# Requires: steps 02 to 12
# ============================================================

source(file.path("01_setup", "01_config.R"))
source(file.path("01_setup", "03_lib_shared.R"))


# Part 1 classifies PLHIV by severity category for 12 indicators under four
# classification scenarios. Part 2 writes the figures and tables. The parts
# run in order.

# ############################################################
# Part 1: single-hazard PLHIV exposure
# 12 Metrics, 4 classification scenarios each
# ############################################################

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(openxlsx)
})

# ============================================================
# 1.0) CONFIGURATION
# ============================================================
config <- list(
  plhiv_files = list(
    ALL = file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", TARGET_YEAR))
  ),

  ewe_metrics = list(
    # --- Lower thresholds (severity) ---
    "TX90p_HWMF" = list(
      path = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWMF/TX90p_ANN_HWMF_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Heatwave"),
    "SPI3_Severity_leq1p0" = list(
      path = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
      flip = TRUE, hazard = "Drought"),
    "R95pTOT" = list(
      path = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_TOT/R95pTOT_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Rainfall"),
    # --- Higher thresholds (severity) ---
    "TX95p_HWMF" = list(
      path = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWMF/TX95p_ANN_HWMF_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Heatwave"),
    "SPI3_Severity_leq1p5" = list(
      path = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR)),
      flip = TRUE, hazard = "Drought"),
    "R99pTOT" = list(
      path = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_TOT/R99pTOT_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Rainfall"),
    # --- Alternative drought indices ---
    "SPEI3_Severity_leq1p0" = list(
      path = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
      flip = TRUE, hazard = "Drought"),
    "SMA3_Severity_leq1p0" = list(
      path = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
      flip = TRUE, hazard = "Drought"),
    "SRI3_Severity_leq1p0" = list(
      path = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
      flip = TRUE, hazard = "Drought"),
    # --- Count variables (lower thresholds) ---
    "TX90p_HWF" = list(
      path = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWF/TX90p_ANN_HWF_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Heatwave"),
    "SPI3_Count_leq1p0" = list(
      path = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Drought"),
    "R95pDAY" = list(
      path = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_DAY/R95pDAY_%d.tif", TARGET_YEAR)),
      flip = FALSE, hazard = "Rainfall")
  ),

  out_dir   = DIR_EXPOSURE,
  # Workbooks and CSV are written to results/exposure/tables and images to
  # results/exposure/figures.
  table_dir = DIR_EXP_TABLES
)

dir.create(config$out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(config$table_dir, showWarnings = FALSE, recursive = TRUE)

# --- fail fast: confirm all EWE rasters exist under results/indices before computing ---
.ewe_paths   <- vapply(config$ewe_metrics, function(m) m$path, character(1))
.ewe_missing <- .ewe_paths[!file.exists(.ewe_paths)]
if (length(.ewe_missing)) {
  message("Missing indicator rasters, all of which step 14 requires:")
  for (p in .ewe_missing) message("   ", p)
  stop("Regenerate/copy the missing EWE layers into results/indices before running step 14.")
}

message("Part 1: single-hazard PLHIV exposure, 12 indicators")


make_break_label <- function(br, discrete = FALSE) {
  k <- length(br) - 1; labs <- character(k)
  for (i in seq_len(k)) {
    if (discrete) {
      lo <- if (i == 1) ceiling(br[1]) else as.integer(ceiling(br[i]))
      hi <- if (i < k) { hi_raw <- br[i+1]; if (hi_raw == floor(hi_raw)) as.integer(hi_raw)-1L else floor(hi_raw) } else floor(br[i+1])
      if (is.na(lo) || is.na(hi)) { labs[i] <- paste0("class ", i); next }
      if (lo > hi) hi <- lo
      labs[i] <- if (lo == hi) as.character(lo) else paste0(lo, "-", hi)
    } else { labs[i] <- paste0(fmt3(br[i]), " to ", fmt3(br[i+1])) }
  }
  labs
}

# ============================================================
# 1.3) CLASSIFICATION SCENARIOS
# ============================================================

# Class boundaries are the median and the two quartiles of the non-zero values
# across the whole domain, so the four classes are quarters of the distribution
# of cells that experienced the event. Zero and missing values are excluded from
# the boundaries and classified as unexposed, because including them would put
# the median at zero for most indicators.
classify_iqr_ssa <- function(v, inside) {
  out <- rep(NA_character_, length(v))
  out[inside & (is.na(v) | v == 0)] <- "Unexposed"
  idx <- which(inside & !is.na(v) & v != 0)
  if (length(idx) == 0) return(list(cats = out, breaks = NULL, labels = NULL))
  vp <- v[idx]; br <- compute_iqr_breaks(vp)
  if (is.null(br) || length(br) < 2) { out[idx] <- "C1"; return(list(cats = out, breaks = br, labels = "all")) }
  cls <- classify_vals(vp, br); out[idx] <- paste0("C", cls)
  list(cats = out, breaks = br, labels = make_break_label(br, is_discrete(vp)))
}

# ============================================================
# 1.4) EXPOSURE CALCULATION
# ============================================================
calculate_exposure <- function(plhiv_draws, cats, r_country_vals, country_levels) {
  PM <- values(plhiv_draws, mat = TRUE); PM[!is.finite(PM)] <- 0; nd <- ncol(PM)
  valid <- !is.na(cats); PMv <- PM[valid, , drop = FALSE]
  cats_v <- cats[valid]; cty_v <- r_country_vals[valid]

  make_summary <- function(idx_local, region_name) {
    if (length(idx_local) == 0) return(NULL)
    cats_sub <- cats_v[idx_local]; PM_sub <- PMv[idx_local, , drop = FALSE]
    plhiv_total <- colSums(PM_sub)
    exp_idx <- which(cats_sub != "Unexposed")
    plhiv_exposed <- if (length(exp_idx) > 0) colSums(PM_sub[exp_idx, , drop = FALSE]) else rep(0, nd)
    all_cats <- c("Unexposed", paste0("C", 1:4), "Any_Exposed")
    results <- list()
    for (cat in all_cats) {
      cat_idx <- if (cat == "Any_Exposed") exp_idx else which(cats_sub == cat)
      if (length(cat_idx) == 0) next
      plhiv_cat <- colSums(PM_sub[cat_idx, , drop = FALSE])
      prop_total <- plhiv_cat / pmax(plhiv_total, 1e-10)
      prop_exposed <- if (cat != "Unexposed" && cat != "Any_Exposed") plhiv_cat / pmax(plhiv_exposed, 1e-10) else rep(NA_real_, nd)
      results[[cat]] <- data.frame(
        region = region_name, category = cat, n_cells = length(cat_idx),
        plhiv_mean = mean(plhiv_cat), plhiv_lower = quantile(plhiv_cat, 0.025),
        plhiv_upper = quantile(plhiv_cat, 0.975),
        pct_of_total_mean = mean(100*prop_total),
        pct_of_total_lower = quantile(100*prop_total, 0.025),
        pct_of_total_upper = quantile(100*prop_total, 0.975),
        pct_of_exposed_mean = mean(100*prop_exposed, na.rm = TRUE),
        pct_of_exposed_lower = quantile(100*prop_exposed, 0.025, na.rm = TRUE),
        pct_of_exposed_upper = quantile(100*prop_exposed, 0.975, na.rm = TRUE),
        stringsAsFactors = FALSE, row.names = NULL)
    }
    bind_rows(results)
  }

  result_ssa <- make_summary(seq_len(nrow(PMv)), "SSA")
  ucty <- sort(unique(cty_v[is.finite(cty_v) & cty_v > 0]))
  result_countries <- bind_rows(lapply(ucty, function(cid) {
    make_summary(which(cty_v == cid), country_name(cid, country_levels))
  }))
  bind_rows(result_ssa, result_countries)
}

# ============================================================
# 1.5) CLASSIFY EVERY INDICATOR
# ============================================================
run_analysis <- function(plhiv_file, cfg) {

  stk <- rast(plhiv_file); tpl <- stk[[1]]
  message(sprintf("  Stack: %d layers, %d cells", nlyr(stk), ncell(tpl)))

  # The analysis domain comes from the draws themselves. A cell is NA in the
  # draws where IHME did not model it.
  plhiv_idx <- grep("^PLHIV_d\\d+$", names(stk))
  inside <- !is.na(values(stk[[plhiv_idx[1]]]))
  message(sprintf("  Mask (from PLHIV draws): %d cells with data, %d NA",
                  sum(inside), sum(!inside)))

  r_country <- find_country_layer(stk)
  r_country_res <- resample(r_country, tpl, method = "near")
  cvals_tpl <- values(r_country_res)
  country_levels <- get_country_levels(r_country_res, r_country)
  message(sprintf("  Countries: %d", nrow(country_levels)))

  plhiv_draws <- stk[[plhiv_idx]]; nd <- nlyr(plhiv_draws)
  message(sprintf("  PLHIV draws: %d", nd))

  all_results <- list(); all_thresholds <- list()

  for (metric_name in names(cfg$ewe_metrics)) {
    mi <- cfg$ewe_metrics[[metric_name]]
    message(sprintf("\n--- %s (%s) ---", metric_name, mi$hazard))
    r_ewe <- rast(mi$path); if (mi$flip) r_ewe <- r_ewe * -1
    r_ewe_crop <- crop(r_ewe, ext(tpl))
    if (!compareGeom(r_ewe_crop, tpl, stopOnError = FALSE))
      r_ewe_crop <- resample(r_ewe_crop, tpl, method = "near")
    # Mask EWE to same domain as PLHIV (using inside vector)
    ewe_vals <- values(r_ewe_crop)
    ewe_vals[!inside] <- NA
    ewe_vals[inside & is.na(ewe_vals)] <- 0
    v <- ewe_vals
    message(sprintf("  Exposed: %d | Unexposed: %d", sum(v[inside] != 0), sum(v[inside] == 0)))

    cls <- classify_iqr_ssa(v, inside)
    exp_one <- calculate_exposure(plhiv_draws, cls$cats, cvals_tpl, country_levels)
    exp_one$metric <- metric_name; exp_one$hazard <- mi$hazard
    all_results <- append(all_results, list(exp_one))
    if (!is.null(cls$breaks))
      all_thresholds <- append(all_thresholds, list(data.frame(
        metric = metric_name, country = "SSA",
        breaks = paste(round(cls$breaks, 4), collapse = " | "),
        labels = paste(cls$labels, collapse = " | "), stringsAsFactors = FALSE)))
  }

  list(results = bind_rows(all_results), thresholds = bind_rows(all_thresholds))
}

# ============================================================
# 1.6) EXCEL OUTPUT
# ============================================================
write_excel <- function(results_df, thresholds_df, out_dir) {
  out_file <- file.path(out_dir, "EWE_exposure_by_severity_class.xlsx")
  wb <- createWorkbook()
  cat_order <- c("Unexposed", "C1", "C2", "C3", "C4", "Any_Exposed")

  add <- function(name, df) {
    addWorksheet(wb, name); writeData(wb, name, df)
    freezePane(wb, name, firstRow = TRUE)
    setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
  }
  add("Exposure", results_df %>%
        arrange(metric, region, factor(category, levels = cat_order)))
  add("Thresholds", thresholds_df)

  saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Written: ", basename(out_file))
}

# ============================================================
# 1.7) RUN PART 1
# ============================================================
result <- run_analysis(config$plhiv_files$ALL, config)
write_excel(result$results, result$thresholds, config$table_dir)

message("Part 1 complete: 12 indicators x 4 classification scenarios")

# ############################################################
# Part 2: figures and tables
# ############################################################

message("\n\n", strrep("=", 80))
message("Part 3: figures and tables")

library(ggplot2)
library(tidyr)
library(scales)
library(sf)

data_dir <- config$table_dir
fig_dir  <- DIR_EXP_FIGS
tbl_dir  <- DIR_EXP_TABLES
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
gadm_path <- FILE_GADM

cat_recode <- c("C1"="Low", "C2"="Moderate", "C3"="High", "C4"="Very High")
cat_levels <- c("Low", "Moderate", "High", "Very High")
class_cols <- c("Low"="#fee8c8", "Moderate"="#fdbb84", "High"="#e34a33", "Very High"="#b30000")

# ============================================================
# 2.1) Load the Part 1 workbook
# ============================================================
.need <- file.path(data_dir, "EWE_exposure_by_severity_class.xlsx")
if (!file.exists(.need))
  stop(basename(.need), " is not in ", data_dir,
       ". The figures are drawn from it, so run this script from the top.")

df <- read.xlsx(.need, sheet = "Exposure")

countries_sf <- st_read(gadm_path, quiet = TRUE)
if (!"COUNTRY" %in% names(countries_sf)) {
  cn_col <- intersect(names(countries_sf), c("COUNTRY","NAME_0","name"))[1]
  if (!is.na(cn_col)) countries_sf$COUNTRY <- countries_sf[[cn_col]]
}

# ============================================================
# 2.2) Figure groups
# ============================================================
# Countries are ordered by the number of PLHIV exposed to the reference indicator.

all_countries <- df %>%
  filter(metric == "TX90p_HWMF", category == "Any_Exposed", region != "SSA") %>%
  arrange(desc(plhiv_mean)) %>% pull(region)
top10 <- head(all_countries, 10)

groups <- list(
  main = list(
    metrics = c("TX90p_HWMF", "R95pTOT", "SPI3_Severity_leq1p0"),
    labels  = c("TX90p_HWMF"           = "a) Heatwave (90th percentile)",
                "R95pTOT"              = "b) Extreme Rainfall (95th percentile)",
                "SPI3_Severity_leq1p0" = "c) Drought (SPI, \u2264-1.0)"),
    slug = "main", title = "Main indicators"),
  sensitivity = list(
    metrics = c("TX95p_HWMF", "R99pTOT", "SPI3_Severity_leq1p5",
                "SPEI3_Severity_leq1p0", "SMA3_Severity_leq1p0", "SRI3_Severity_leq1p0"),
    labels  = c("TX95p_HWMF"            = "a) Heatwave\n(95th percentile)",
                "R99pTOT"               = "b) Extreme Rainfall\n(99th percentile)",
                "SPI3_Severity_leq1p5"  = "c) Drought\n(SPI, \u2264-1.5)",
                "SPEI3_Severity_leq1p0" = "d) Drought\n(SPEI, \u2264-1.0)",
                "SMA3_Severity_leq1p0"  = "e) Drought\n(SMA, \u2264-1.0)",
                "SRI3_Severity_leq1p0"  = "f) Drought\n(SRI, \u2264-1.0)"),
    slug = "sensitivity", title = "Sensitivity indicators")
)

# ============================================================
# 2.3) Figures
# ============================================================
# Share of each country's PLHIV in each severity class, class boundaries taken
# across the whole of SSA. One panel per indicator, all 43 countries.
fig_severity_by_country <- function(g) {
  fn <- paste0("exposure_severity_by_country_", g$slug)
  d <- df %>% filter(metric %in% g$metrics, region %in% all_countries,
                     category %in% paste0("C", 1:4)) %>%
    mutate(Panel = factor(g$labels[metric], levels = g$labels),
           Severity = factor(cat_recode[category], levels = cat_levels),
           region = factor(region, levels = rev(all_countries)))
  if (nrow(d) == 0) { message("Skipped: ", fn); return(invisible()) }
  totals <- df %>% filter(metric %in% g$metrics, region %in% all_countries,
                          category == "Any_Exposed") %>%
    mutate(Panel = factor(g$labels[metric], levels = g$labels),
           region = factor(region, levels = rev(all_countries)),
           ann = paste0(fmt_num(plhiv_mean), " (", sprintf("%.0f%%", pct_of_total_mean), ")"))
  nc <- length(all_countries)
  p <- ggplot(d, aes(x = pct_of_total_mean, y = region, fill = Severity)) +
    geom_col() +
    geom_text(data = totals, aes(x = pct_of_total_mean, y = region, label = ann, fill = NULL),
              hjust = -0.05, size = 1.6, colour = "grey30") +
    scale_fill_manual(values = class_cols, name = "Severity") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35)), labels = label_percent(scale = 1)) +
    labs(x = NULL, y = "% of national PLHIV") +
    facet_wrap(~ Panel, scales = "free_x", ncol = length(g$metrics)) +
    theme_ewe(base_size = 9) +
    theme(panel.spacing = unit(0.8, "lines"), axis.text.y = element_text(size = 6))
  ggsave(file.path(fig_dir, paste0(fn, ".jpeg")), p,
         width = max(14, 3 + 3.5 * length(g$metrics)), height = max(6, nc * 0.25 + 2), dpi = 600)
  message(fn)
}

# Share of each country's PLHIV falling in the Very High class, as a map.
fig_very_high_map <- function(g) {
  fn <- paste0("exposure_very_high_map_", g$slug)
  d <- df %>% filter(metric %in% g$metrics, category == "C4", region != "SSA") %>%
    mutate(Panel = factor(g$labels[metric], levels = g$labels))
  if (nrow(d) == 0) { message("Skipped: ", fn); return(invisible()) }
  map_sf <- countries_sf %>% left_join(d, by = c("COUNTRY" = "region"))
  drawable <- drawable_or_skip(map_sf, "Panel", fn)
  if (is.null(drawable)) return(invisible())
  p <- ggplot(drawable) +
    geom_sf(data = countries_sf, fill = "grey90", colour = "grey60", linewidth = 0.15) +
    geom_sf(aes(fill = pct_of_total_mean), colour = "grey40", linewidth = 0.2) +
    scale_fill_distiller(palette = "Reds", direction = 1, na.value = "grey90",
                         labels = label_percent(scale = 1), name = "% PLHIV in\nVery High") +
    facet_wrap(~ Panel, ncol = length(g$metrics)) +
    theme_map()
  ggsave(file.path(fig_dir, paste0(fn, ".jpeg")), p,
         width = max(14, 4 + 3.5 * length(g$metrics)), height = 5.5, dpi = 600)
  message(fn)
}

for (g in groups) { fig_severity_by_country(g); fig_very_high_map(g) }

# ============================================================
# 2.4) Tables
# ============================================================
all_metrics <- unlist(lapply(groups, `[[`, "metrics"), use.names = FALSE)
tbl_lab <- gsub("\n", " ", unlist(lapply(groups, `[[`, "labels")))

tbl1 <- df %>%
  filter(metric %in% all_metrics, region == "SSA",
         category %in% c("Any_Exposed", paste0("C", 1:4))) %>%
  mutate(Hazard = tbl_lab[metric], Metric = metric,
         Category = ifelse(category == "Any_Exposed", "Any exposure", cat_recode[category]),
         Category = factor(Category, levels = c("Any exposure", cat_levels)),
         Cells = format(n_cells, big.mark = ","),
         `PLHIV (95% UI)` = sprintf("%.2fM [%.2fM-%.2fM]",
                                    plhiv_mean/1e6, plhiv_lower/1e6, plhiv_upper/1e6),
         `% of total PLHIV` = sprintf("%.1f%% [%.1f%%-%.1f%%]",
                                      pct_of_total_mean, pct_of_total_lower, pct_of_total_upper),
         `% of exposed PLHIV` = ifelse(is.na(pct_of_exposed_mean), "--",
                                       sprintf("%.1f%% [%.1f%%-%.1f%%]", pct_of_exposed_mean,
                                               pct_of_exposed_lower, pct_of_exposed_upper))) %>%
  select(Hazard, Metric, Category, Cells, `PLHIV (95% UI)`,
         `% of total PLHIV`, `% of exposed PLHIV`) %>%
  arrange(Metric, Category)
write.csv(tbl1, file.path(tbl_dir, "exposure_by_severity_ssa.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

tbl2 <- df %>%
  filter(metric %in% groups$main$metrics, region %in% top10,
         category %in% c("Any_Exposed", "C4")) %>%
  mutate(Hazard = tbl_lab[metric], Country = region,
         Category = ifelse(category == "Any_Exposed", "Any exposure", "Very High"),
         `PLHIV (95% UI)` = sprintf("%.0f [%.0f-%.0f]", plhiv_mean, plhiv_lower, plhiv_upper),
         `% of national PLHIV` = sprintf("%.1f%%", pct_of_total_mean)) %>%
  select(Hazard, Country, Category, n_cells, `PLHIV (95% UI)`, `% of national PLHIV`) %>%
  arrange(Hazard, desc(Category == "Any exposure"), Country)
write.csv(tbl2, file.path(tbl_dir, "exposure_top10_countries.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("Done. Next: step 15.")
