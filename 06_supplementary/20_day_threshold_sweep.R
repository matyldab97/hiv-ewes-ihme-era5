# ============================================================
# Supplementary: exposure against the day threshold
# ------------------------------------------------------------
# Step 15 reports exposure at a single cut-off of seven days. This step
# recomputes it at every cut-off from one to fourteen days, so that the reported
# figure can be read against the whole curve.
#
# Only the two day-based EWEs are checked, heatwave days and extreme rainfall
# days. Drought is counted in months, so a day threshold does not apply to it.
#
# Summarised over the 1000 PLHIV draws, for sub-Saharan Africa and per country,
# reporting the mean and a 95% interval.
#
# Requires: steps 05, 12 and the heatwave step 04.
# Output: results/supplementary/day_threshold_sweep
# ============================================================

source(file.path("01_setup", "01_config.R"))
source(file.path("01_setup", "03_lib_shared.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(openxlsx); library(ggplot2)
  library(sf); library(scales)
})

out_dir <- file.path(DIR_SUPP, "day_threshold_sweep")
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

yr         <- TARGET_YEAR
THRESHOLDS <- 1:14

plhiv_file <- file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", yr))
hazards <- list(
  HW = list(path = file.path(DIR_INDICES, "Heatwaves", "TX90p", "ANN", "HWF",
                             sprintf("TX90p_ANN_HWF_%d.tif", yr)),
            label = "Heatwave (90th percentile)"),
  RF = list(path = file.path(DIR_INDICES, "R95p", "ANNUAL_DAY",
                             sprintf("R95pDAY_%d.tif", yr)),
            label = "Extreme rainfall (95th percentile)")
)

for (f in c(plhiv_file, vapply(hazards, `[[`, character(1), "path")))
  if (!file.exists(f)) stop("Missing input, run the earlier steps first:\n  ", f)

message("Loading PLHIV draws ...")
stk <- rast(plhiv_file)
plhiv_draws <- stk[[grep("^PLHIV_d[0-9]+$", names(stk))]]
if (nlyr(plhiv_draws) == 0) stop("No PLHIV_dN layers in ", plhiv_file)
nd     <- nlyr(plhiv_draws)
tpl    <- plhiv_draws[[1]]
inside <- !is.na(values(tpl))
message(sprintf("  %d cells, %d draws, %d cells with data", ncell(tpl), nd, sum(inside)))

r_country <- find_country_layer(stk)
if (!compareGeom(r_country, tpl, stopOnError = FALSE))
  r_country <- resample(r_country, tpl, method = "near")
country_levels <- get_country_levels(r_country)

PM <- values(plhiv_draws, mat = TRUE); PM[!is.finite(PM)] <- 0
PMv   <- PM[inside, , drop = FALSE]
cty_v <- values(r_country)[inside]
ucty  <- sort(unique(cty_v[is.finite(cty_v) & cty_v > 0]))

# The denominator does not depend on the threshold, so the per-region totals
# are computed once rather than 14 times per hazard.
region_idx <- c(list(SSA = seq_len(nrow(PMv))),
                setNames(lapply(ucty, function(cid) which(cty_v == cid)),
                         vapply(ucty, function(cid) country_name(cid, country_levels),
                                character(1))))
region_tot <- lapply(region_idx, function(ix) colSums(PMv[ix, , drop = FALSE]))

align_days <- function(path) {
  r <- crop(rast(path), ext(tpl))
  if (!compareGeom(r, tpl, stopOnError = FALSE)) r <- resample(r, tpl, method = "near")
  v <- values(r); v[!inside] <- NA; v[inside & is.na(v)] <- 0
  v[inside]
}

rows <- list()
for (hz in names(hazards)) {
  h  <- hazards[[hz]]
  dv <- align_days(h$path)
  message(sprintf("\n%s: sweeping 1 to %d days", h$label, max(THRESHOLDS)))
  for (N in THRESHOLDS) {
    exp_v <- dv >= N
    for (rg in names(region_idx)) {
      ix <- region_idx[[rg]]; es <- exp_v[ix]
      pe <- if (any(es)) colSums(PMv[ix, , drop = FALSE][es, , drop = FALSE]) else rep(0, nd)
      prop <- 100 * pe / pmax(region_tot[[rg]], 1e-10)
      rows[[length(rows) + 1]] <- data.frame(
        hazard = h$label, threshold_days = N, region = rg, n_cells_exposed = sum(es),
        plhiv_exposed_mean   = mean(pe),
        plhiv_exposed_lower  = quantile(pe, 0.025),
        plhiv_exposed_upper  = quantile(pe, 0.975),
        pct_of_total_mean    = mean(prop),
        pct_of_total_lower   = quantile(prop, 0.025),
        pct_of_total_upper   = quantile(prop, 0.975),
        stringsAsFactors = FALSE, row.names = NULL)
    }
    message(sprintf("  >= %2d days", N))
  }
}
sweep <- bind_rows(rows)

ssa_tab <- sweep %>% filter(region == "SSA") %>% arrange(hazard, threshold_days)
cty_tab <- sweep %>% filter(region != "SSA") %>% arrange(hazard, region, threshold_days)
write.csv(ssa_tab, file.path(out_dir, sprintf("day_threshold_1to%d_SSA.csv", max(THRESHOLDS))),
          row.names = FALSE)
write.csv(cty_tab, file.path(out_dir, sprintf("day_threshold_1to%d_Country.csv", max(THRESHOLDS))),
          row.names = FALSE)
wb <- createWorkbook()
addWorksheet(wb, "SSA");     writeData(wb, "SSA", ssa_tab);     freezePane(wb, "SSA", firstRow = TRUE)
addWorksheet(wb, "Country"); writeData(wb, "Country", cty_tab); freezePane(wb, "Country", firstRow = TRUE)
saveWorkbook(wb, file.path(out_dir, sprintf("day_threshold_1to%d_exposure.xlsx", max(THRESHOLDS))),
             overwrite = TRUE)
message("\nTables written to ", out_dir)

# The value step 15 reports, marked on the curve so the two cannot drift apart.
p_ssa <- ggplot(ssa_tab, aes(threshold_days, pct_of_total_mean, colour = hazard, fill = hazard)) +
  geom_ribbon(aes(ymin = pct_of_total_lower, ymax = pct_of_total_upper),
              alpha = 0.15, colour = NA) +
  geom_vline(xintercept = 7, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = THRESHOLDS) +
  labs(x = "Exposure threshold (at least N days)", y = "% of PLHIV exposed",
       colour = NULL, fill = NULL,
       title = "Exposed proportion against the day threshold",
       subtitle = "Dashed line marks the seven-day cut-off reported in step 15") +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(out_dir, "figures", "day_threshold_SSA_trend.jpeg"), p_ssa,
       width = 9, height = 6, dpi = 300)

for (hz in names(hazards)) {
  h  <- hazards[[hz]]
  dd <- cty_tab %>% filter(hazard == h$label)
  p <- ggplot(dd, aes(threshold_days, pct_of_total_mean)) +
    geom_ribbon(aes(ymin = pct_of_total_lower, ymax = pct_of_total_upper),
                alpha = 0.2, fill = "steelblue") +
    geom_line(colour = "steelblue", linewidth = 0.6) +
    geom_vline(xintercept = 7, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    facet_wrap(~ region, ncol = 6) +
    scale_x_continuous(breaks = c(1, 7, 14)) +
    labs(x = "at least N days", y = "% of national PLHIV exposed",
         title = paste0(h$label, ": exposed proportion against the day threshold, by country")) +
    theme_bw(base_size = 8)
  ggsave(file.path(out_dir, "figures", paste0("day_threshold_country_trend_", hz, ".jpeg")),
         p, width = 14, height = 10, dpi = 300)
}

countries_sf <- st_read(FILE_GADM, quiet = TRUE)
if (!"COUNTRY" %in% names(countries_sf)) {
  cn <- intersect(names(countries_sf), c("COUNTRY","NAME_0","name"))[1]
  if (!is.na(cn)) countries_sf$COUNTRY <- countries_sf[[cn]]
}

for (hz in names(hazards)) {
  h  <- hazards[[hz]]
  lab <- function(n) paste0("at least ", n, " day", ifelse(n == 1, "", "s"))
  dd <- cty_tab %>% filter(hazard == h$label) %>%
    mutate(thr_lab = factor(lab(threshold_days), levels = lab(THRESHOLDS)))
  map_sf <- countries_sf %>% left_join(dd, by = c("COUNTRY" = "region"))
  drawable <- drawable_or_skip(map_sf, "thr_lab", paste0("day_threshold_MAPS_", hz))
  if (is.null(drawable)) next
  p <- ggplot() +
    geom_sf(data = countries_sf, fill = "grey90", colour = "grey60", linewidth = 0.08) +
    geom_sf(data = drawable, aes(fill = pct_of_total_mean), colour = "grey40", linewidth = 0.08) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, na.value = "grey90",
                         limits = c(0, 100), labels = label_percent(scale = 1),
                         name = "% national\nPLHIV exposed") +
    facet_wrap(~ thr_lab, ncol = 5) +
    labs(title = paste0(h$label, ": proportion of PLHIV exposed, by day threshold")) +
    theme_map()
  ggsave(file.path(out_dir, "figures", paste0("day_threshold_MAPS_", hz, ".jpeg")),
         p, width = 16, height = 11, dpi = 300)
  message("  day_threshold_MAPS_", hz)
}

message("\nDone. Output: ", out_dir)
