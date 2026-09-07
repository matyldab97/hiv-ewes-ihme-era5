# ============================================================
# Threshold exposure
# ------------------------------------------------------------
# PLHIV living in cells with at least seven heatwave days, at least seven
# extreme rainfall days, or at least one drought month.
#
# Requires: steps 02 to 12
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))
source(file.path("R", "01_setup", "03_lib_shared.R"))


# Reads the PLHIV draws with the country layer, the heatwave-day and
# extreme-rainfall-day surfaces for the target year, and the GADM polygons.
#
# Output: results/exposure, one workbook with a continent-wide sheet and a
# per-country sheet, the corresponding CSV tables, and the figures.

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(openxlsx)
  library(ggplot2); library(scales); library(sf)
})

# ============================================================
# 0) CONFIG - edit paths here only
# ============================================================
cfg <- list(

  # The Monte Carlo draw stack: PLHIV_d1, PLHIV_d2, ... plus a categorical
  # 'country' layer, as step 12 writes it.
  plhiv_files = list(
    ALL = file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", TARGET_YEAR))
  ),

  # Hazard rasters for the target year. Each hazard carries its own
  # threshold: heatwave and rainfall use >= 7 days, drought uses >= 1 month.
  hazards = list(
    HW = list(
      path      = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWF/TX90p_ANN_HWF_%d.tif", TARGET_YEAR)),
      hazard    = "Heatwave",
      panel     = "a) Heatwave (90th percentile)\n\u22657 days",
      threshold = 7, unit = "days"
    ),
    RF = list(
      path      = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_DAY/R95pDAY_%d.tif", TARGET_YEAR)),
      hazard    = "Extreme rainfall",
      panel     = "b) Extreme Rainfall (95th percentile)\n\u22657 days",
      threshold = 7, unit = "days"
    ),
    DR = list(
      # Number of months in the target year with SPI-3 <= -1.0.
      path      = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
      hazard    = "Drought",
      panel     = "c) Drought (SPI, \u2264-1.0)\n\u22651 month",
      threshold = 1, unit = "months"
    ),

    ## ---- 6 SENSITIVITY indicators (for supplement S3/S10) ----
    HW95 = list(
      path      = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWF/TX95p_ANN_HWF_%d.tif", TARGET_YEAR)),
      hazard    = "Heatwave (95th percentile)",
      panel     = "a) Heatwave (95th percentile)\n\u22657 days",
      threshold = 7, unit = "days", sens = TRUE
    ),
    RF99 = list(
      path      = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_DAY/R99pDAY_%d.tif", TARGET_YEAR)),
      hazard    = "Extreme rainfall (99th percentile)",
      panel     = "b) Extreme Rainfall (99th percentile)\n\u22657 days",
      threshold = 7, unit = "days", sens = TRUE
    ),
    DR_SPI15 = list(
      path      = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_month_count_leq-1.5.tif", TARGET_YEAR)),
      hazard    = "Drought (SPI, <=-1.5)",
      panel     = "c) Drought (SPI, \u2264-1.5)\n\u22651 month",
      threshold = 1, unit = "months", sens = TRUE
    ),
    DR_SPEI = list(
      path      = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
      hazard    = "Drought (SPEI)",
      panel     = "d) Drought (SPEI, \u2264-1.0)\n\u22651 month",
      threshold = 1, unit = "months", sens = TRUE
    ),
    DR_SRI = list(
      path      = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
      hazard    = "Drought (SRI)",
      panel     = "e) Drought (SRI, \u2264-1.0)\n\u22651 month",
      threshold = 1, unit = "months", sens = TRUE
    ),
    DR_SMA = list(
      path      = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_month_count_leq-1.0.tif", TARGET_YEAR)),
      hazard    = "Drought (SMA)",
      panel     = "f) Drought (SMA, \u2264-1.0)\n\u22651 month",
      threshold = 1, unit = "months", sens = TRUE
    )
  ),

  # Country polygons (for maps)
  gadm_path = FILE_GADM,

  # Output directory 
  out_dir = DIR_EXPOSURE
)

fig_dir <- DIR_EXP_FIGS      # images only
tbl_dir <- DIR_EXP_TABLES    # workbooks and CSV only
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir,     showWarnings = FALSE, recursive = TRUE)

message("PLHIV EXPOSED AT OR ABOVE A PER-HAZARD THRESHOLD")
for (.hz in cfg$hazards)
  message(sprintf("  %-34s >= %d %s%s", .hz$hazard, .hz$threshold, .hz$unit,
                  if (isTRUE(.hz$sens)) "   (sensitivity)" else ""))
message(sprintf("Output: %s", cfg$out_dir))

# ============================================================
# 1) EXPOSURE AT THE PER-HAZARD THRESHOLD
# ============================================================
run_threshold_exposure <- function(plhiv_file, cfg) {
  stk <- rast(plhiv_file); tpl <- stk[[1]]

  # Mask straight from PLHIV draws
  plhiv_idx <- grep("^PLHIV_d\\d+$", names(stk))
  if (length(plhiv_idx) == 0)
    stop(sprintf("No PLHIV_dN layers found in %s", plhiv_file))
  inside <- !is.na(values(stk[[plhiv_idx[1]]]))
  message(sprintf("  Mask: %d cells with PLHIV data", sum(inside)))

  # Country layer
  r_country     <- find_country_layer(stk)
  r_country_res <- resample(r_country, tpl, method = "near")
  cvals         <- values(r_country_res)
  country_levels <- get_country_levels(r_country_res, r_country)
  message(sprintf("  Countries detected: %d", nrow(country_levels)))

  # PLHIV draws matrix
  plhiv_draws <- stk[[plhiv_idx]]
  nd <- nlyr(plhiv_draws)
  PM <- values(plhiv_draws, mat = TRUE)
  PM[!is.finite(PM)] <- 0
  PMv   <- PM[inside, , drop = FALSE]
  cty_v <- cvals[inside]
  message(sprintf("  Draws: %d", nd))

  # Per-hazard binary >= threshold
  summarise_hazard <- function(hz_key) {
    hz <- cfg$hazards[[hz_key]]
    r  <- rast(hz$path)
    rc <- crop(r, ext(tpl))
    if (!compareGeom(rc, tpl, stopOnError = FALSE))
      rc <- resample(rc, tpl, method = "near")
    v <- values(rc)
    v[!inside] <- NA
    v[inside & is.na(v)] <- 0
    exposed_cell <- inside & v >= hz$threshold
    exp_v <- exposed_cell[inside]
    message(sprintf("  %-34s >= %d %-6s: %d exposed cells (of %d)",
                    hz$hazard, hz$threshold, hz$unit, sum(exp_v), sum(inside)))

    one_region <- function(idx_local, region_name) {
      if (length(idx_local) == 0) return(NULL)
      exp_sub   <- exp_v[idx_local]
      PM_sub    <- PMv[idx_local, , drop = FALSE]
      plhiv_tot <- colSums(PM_sub)
      plhiv_exp <- if (sum(exp_sub) > 0) colSums(PM_sub[exp_sub, , drop = FALSE]) else rep(0, nd)
      prop      <- plhiv_exp / pmax(plhiv_tot, 1e-10)
      data.frame(
        hazard               = hz$hazard,
        panel                = hz$panel,
        sens                 = isTRUE(hz$sens),
        threshold            = hz$threshold,
        threshold_unit       = hz$unit,
        region               = region_name,
        n_cells_total        = length(idx_local),
        n_cells_exposed      = sum(exp_sub),
        plhiv_exposed_mean = mean(plhiv_exp),
        plhiv_exposed_lower  = quantile(plhiv_exp, 0.025),
        plhiv_exposed_upper  = quantile(plhiv_exp, 0.975),
        plhiv_total_mean   = mean(plhiv_tot),
        plhiv_total_lower    = quantile(plhiv_tot, 0.025),
        plhiv_total_upper    = quantile(plhiv_tot, 0.975),
        pct_of_total_mean  = mean(100 * prop),
        pct_of_total_lower   = quantile(100 * prop, 0.025),
        pct_of_total_upper   = quantile(100 * prop, 0.975),
        stringsAsFactors = FALSE, row.names = NULL
      )
    }

    ssa <- one_region(seq_len(nrow(PMv)), "SSA")
    ucty <- sort(unique(cty_v[is.finite(cty_v) & cty_v > 0]))
    cty_res <- bind_rows(lapply(ucty, function(cid) {
      one_region(which(cty_v == cid), country_name(cid, country_levels))
    }))
    bind_rows(ssa, cty_res)
  }

  bind_rows(lapply(names(cfg$hazards), summarise_hazard))
}

# ============================================================
# 2) Run and write the workbook
# ============================================================
all_res <- run_threshold_exposure(cfg$plhiv_files$ALL, cfg)

wb <- createWorkbook()
add_sheet <- function(name, df) {
  addWorksheet(wb, name); writeData(wb, name, df)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("SSA",     all_res %>% filter(region == "SSA") %>% arrange(hazard))
add_sheet("Country", all_res %>% filter(region != "SSA") %>% arrange(hazard, region))
out_file <- file.path(DIR_EXP_TABLES, "threshold_exposure.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Written: ", basename(out_file))

# ============================================================
# 3) Figures
# ============================================================
message("\n--- Figures ---")

df_7d_all <- all_res
df_7d      <- df_7d_all %>% filter(!sens)   # 3 main hazards, manuscript map
df_7d_sens <- df_7d_all %>% filter(sens)    # 6 sensitivity hazards, supplement map
df_7d_ssa  <- df_7d %>% filter(region == "SSA")
df_7d_cty  <- df_7d %>% filter(region != "SSA")

## panel factor levels: main (3) and sensitivity (6), taken in config order
main_panels <- vapply(Filter(function(h) !isTRUE(h$sens), cfg$hazards), `[[`, character(1), "panel", USE.NAMES = FALSE)
sens_panels <- vapply(Filter(function(h)  isTRUE(h$sens), cfg$hazards), `[[`, character(1), "panel", USE.NAMES = FALSE)
panel_levels <- main_panels

countries_sf <- st_read(cfg$gadm_path, quiet = TRUE)
if (!"COUNTRY" %in% names(countries_sf)) {
  cn_col <- intersect(names(countries_sf), c("COUNTRY", "NAME_0", "name"))[1]
  if (!is.na(cn_col)) countries_sf$COUNTRY <- countries_sf[[cn_col]]
}

# ---- Choropleth map, three main hazards -------------------------------------
dM <- df_7d_cty %>% mutate(panel = factor(panel, levels = panel_levels))
map_sf <- countries_sf %>% left_join(dM, by = c("COUNTRY" = "region"))
drawable_m <- drawable_or_skip(map_sf, "panel", "sevenday_exposure_map_main")
if (!is.null(drawable_m)) {
pM <- ggplot(drawable_m) +
  geom_sf(data = countries_sf, fill = "grey90", colour = "grey60", linewidth = 0.15) +
  geom_sf(aes(fill = pct_of_total_mean), colour = "grey40", linewidth = 0.2) +
  facet_wrap(~ panel, ncol = 3) +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "grey90",
                       limits = c(0, 100),
                       labels = label_percent(scale = 1),
                       name = "% PLHIV\nexposed") +
  theme_map()

ggsave(file.path(fig_dir, "sevenday_exposure_map_main.jpeg"), pM,
       width = 15, height = 5.5, dpi = 600, type = "cairo")
ggsave(file.path(fig_dir, "sevenday_exposure_map_main.pdf"),  pM,
       width = 15, height = 5.5, device = cairo_pdf)
message("  sevenday_exposure_map_main")
}

# ---- Choropleth map, six sensitivity hazards (supplement) -------------------
if (nrow(df_7d_sens) > 0) {
  dMs <- df_7d_sens %>% mutate(panel = factor(panel, levels = sens_panels))
  map_sf_s <- countries_sf %>% left_join(dMs, by = c("COUNTRY" = "region"))
  drawable_s <- drawable_or_skip(map_sf_s, "panel", "sevenday_exposure_map_sensitivity")
  if (!is.null(drawable_s)) {
  pMs <- ggplot(drawable_s) +
    geom_sf(data = countries_sf, fill = "grey90", colour = "grey60", linewidth = 0.15) +
    geom_sf(aes(fill = pct_of_total_mean), colour = "grey40", linewidth = 0.2) +
    facet_wrap(~ panel, ncol = 3) +
    scale_fill_distiller(palette = "Reds", direction = 1, na.value = "grey90",
                         limits = c(0, 100), labels = label_percent(scale = 1),
                         name = "% PLHIV\nexposed") +
    theme_map()
  ggsave(file.path(fig_dir, "sevenday_exposure_map_sensitivity.jpeg"), pMs,
         width = 15, height = 10, dpi = 600, type = "cairo")
  ggsave(file.path(fig_dir, "sevenday_exposure_map_sensitivity.pdf"),  pMs,
         width = 15, height = 10, device = cairo_pdf)
  message("  sevenday_exposure_map_sensitivity")
  }
}

# ============================================================
# 4) Tables
# ============================================================
tbl1 <- df_7d_ssa %>%
  mutate(Hazard = panel,
         Cells  = format(n_cells_exposed, big.mark = ","),
         `PLHIV exposed (95% UI)` = sprintf("%.2fM [%.2fM\u2013%.2fM]",
                                            plhiv_exposed_mean / 1e6,
                                            plhiv_exposed_lower  / 1e6,
                                            plhiv_exposed_upper  / 1e6),
         `% of total PLHIV (95% UI)` = sprintf("%.1f%% [%.1f%%\u2013%.1f%%]",
                                                pct_of_total_mean,
                                                pct_of_total_lower,
                                                pct_of_total_upper)) %>%
  select(Hazard, Cells,
         `PLHIV exposed (95% UI)`,
         `% of total PLHIV (95% UI)`)
write.csv(tbl1, file.path(tbl_dir, "threshold_exposure_ssa.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


tbl2 <- df_7d_cty %>%
  arrange(hazard, desc(plhiv_exposed_mean)) %>%
  mutate(Hazard  = panel,
         Country = region,
         `PLHIV exposed (95% UI)` = sprintf("%.0f [%.0f\u2013%.0f]",
                                            plhiv_exposed_mean,
                                            plhiv_exposed_lower,
                                            plhiv_exposed_upper),
         `% of national PLHIV (95% UI)` = sprintf("%.1f%% [%.1f%%\u2013%.1f%%]",
                                                    pct_of_total_mean,
                                                    pct_of_total_lower,
                                                    pct_of_total_upper)) %>%
  select(Hazard, Country, n_cells_exposed,
         `PLHIV exposed (95% UI)`,
         `% of national PLHIV (95% UI)`)
write.csv(tbl2, file.path(tbl_dir, "threshold_exposure_by_country.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("Done. Next: step 16.")
