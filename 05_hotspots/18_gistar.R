# ============================================================
# Getis-Ord Gi*
# ------------------------------------------------------------
# Local Gi* (Getis and Ord, 1992, Geographical Analysis 24:189,
# doi 10.1111/j.1538-4632.1992.tb00261.x) for every input surface,
# continent-wide and within each country, with Benjamini-Hochberg correction
# (1995, JRSS B 57:289) of the local p-values.
#
# Category codes give the strictest level at which a cell survives correction:
# 1 for q <= 0.10, 2 for q <= 0.05, 3 for q <= 0.01, signed by the z-score.
#
# Requires: step 16
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))
source(file.path("R", "01_setup", "03_lib_shared.R"))

# Twelve input surfaces are analysed: the nine indicators and three summaries
# of the prevalence draws. Each is analysed continent-wide and separately
# within each country.
#
# Input:  the transformed stack from step 16
# Output: results/gistar, one raster of z-scores and p-values and one of
#         category codes per surface, the significance and PLHIV-in-hotspot
#         tables, and the maps

suppressPackageStartupMessages({
  library(terra)
  library(spdep)
  library(sf)
})

# ============================================================
# CONFIG
# ============================================================

config <- list(
  in_dir     = file.path(DIR_GI_IN, "ready"),
  in_pattern = "_GiStar_ready\\.tif$",
  out_dir    = DIR_GI_OUT,

  # PLHIV source for counting people in hotspots
  plhiv_ui_path = file.path(DIR_MC, sprintf("ALL_15_59_%d_UI_0p25.tif", TARGET_YEAR)),
  plhiv_draws_path = file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", TARGET_YEAR)),

  out_suffix_results     = "_GiStar_results.tif",
  out_suffix_categorical = "_GiStar_categorical.tif",

  # Gi* rather than Gi, so a cell is part of its own neighbourhood
  neigh_type   = "queen",
  style        = "B",          
  include_self = TRUE,

  # One test per land cell is roughly 45,000 tests per surface, so the
  # uncorrected p-values are not usable on their own. Benjamini-Hochberg
  # controls the false discovery rate.
  fdr_methods = c("BH"),

  # Significance bands, strictly decreasing.
  # Category code k (and -k) == cell survives at gi_alphas[k].
  gi_alphas   = c("90%" = 0.10, "95%" = 0.05, "99%" = 0.01),

  # The level that counts as "significant" for binary sig layers,
  # overlap maps and the hotspot population counts. Must be one of gi_alphas.
  fdr_alpha   = 0.05,

  ewe_layers = c(
    "TX90p_HWMF", "TX95p_HWMF",
    "R95pTOT", "R99pTOT",
    "SPI3_leq1p0", "SPI3_leq1p5", "SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0"
  ),

  prev_draw_pattern = "^PREVpct_d[0-9]+$",
  min_cells_country = 1,

  gadm_path = FILE_GADM,

  country_names_csv = NULL,

  var_labels = c(
    TX90p_HWMF    = "Heatwave (90th percentile)",
    R95pTOT       = "Extreme Rainfall (95th percentile)",
    SPI3_leq1p0   = "Drought (SPI, \u2264-1.0)",
    TX95p_HWMF    = "Heatwave (95th percentile)",
    R99pTOT       = "Extreme Rainfall (99th percentile)",
    SPI3_leq1p5   = "Drought (SPI, \u2264-1.5)",
    SPEI3_leq1p0  = "Drought (SPEI, \u2264-1.0)",
    SRI3_leq1p0   = "Drought (SRI, \u2264-1.0)",
    SMA3_leq1p0   = "Drought (SMA, \u2264-1.0)",
    PREVpct_mean  = "HIV Prevalence (mean)",
    PREVpct_lower = "HIV Prevalence (lower)",
    PREVpct_upper = "HIV Prevalence (upper)"
  ),

  main_vars = c("TX90p_HWMF", "R95pTOT", "SPI3_leq1p0")
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TMP, showWarnings = FALSE)
terraOptions(progress = 1)

# Z cut-offs derived from gi_alphas
Z_CUTS <- qnorm(1 - config$gi_alphas / 2)   # length 3, increasing

stopifnot(all(diff(config$gi_alphas) < 0))

# alpha -> category code magnitude (0.10 -> 1, 0.05 -> 2, 0.01 -> 3)
alpha_to_code <- function(alpha) {
  k <- which(abs(config$gi_alphas - alpha) < 1e-12)
  if (!length(k)) stop("fdr_alpha (", alpha, ") must be one of config$gi_alphas")
  k
}
SIG_CODE <- alpha_to_code(config$fdr_alpha)   # 2 when fdr_alpha = 0.05

# Both tables are rewritten from scratch below. Removing them first means a run
# that fails part way through leaves no table at all, rather than one from a
# previous run that could be mistaken for the current result.
stale <- c("PLHIV_Gi_multithreshold_hotcold.csv",
           "PLHIV_in_hotspots_coldspots.csv")
for (f in file.path(config$out_dir, stale)) {
  if (file.exists(f)) { unlink(f); message("  Removed stale CSV: ", basename(f)) }
}

message("  Getis-Ord Gi* analysis, Benjamini-Hochberg correction")
message("  Surfaces: 12 (9 indicators, 3 prevalence summaries) | continent-wide and per country")

# ============================================================
# Neighbour functions
# ============================================================

build_nb_ssa <- function(r_stack) {
  valid <- !is.na(values(r_stack$country, mat = FALSE))
  nb <- spdep::cell2nb(nrow(r_stack), ncol(r_stack),
                       type = config$neigh_type, torus = FALSE)
  nb <- spdep::subset.nb(nb, valid)
  if (config$include_self) nb <- spdep::include.self(nb)
  message(sprintf("  SSA neighbour graph: %s land cells",
                  format(sum(valid), big.mark = ",")))
  list(nb = nb, valid = valid)
}

build_nb_country <- function(r_stack) {
  message("  Precomputing country neighbour structures...")
  t0 <- Sys.time()
  country_vals <- values(r_stack$country, mat = FALSE)
  valid <- !is.na(country_vals)
  nb_ssa <- spdep::cell2nb(nrow(r_stack), ncol(r_stack),
                           type = config$neigh_type, torus = FALSE)
  nb_ssa <- spdep::subset.nb(nb_ssa, valid)
  if (config$include_self) nb_ssa <- spdep::include.self(nb_ssa)
  country_land <- country_vals[valid]
  sizes <- table(country_land)
  valid_ids <- as.numeric(names(sizes)[sizes >= config$min_cells_country])
  structs <- list()
  for (cid in valid_ids) {
    mask <- (country_land == cid)
    structs[[as.character(cid)]] <- list(
      country_id = cid, mask_land = mask,
      nb = spdep::subset.nb(nb_ssa, mask),
      indices_land = which(mask), n_cells = sum(mask))
  }
  message(sprintf("    %d countries in %.1f sec",
                  length(structs),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(valid = valid, country_vals_land = country_land,
       country_structures = structs)
}

# ============================================================
# Gi* CORE
# ============================================================

gi_ssa <- function(r_layer, nb, valid) {
  x <- values(r_layer, mat = FALSE)
  keep <- is.finite(x[valid])
  out <- rast(r_layer)
  if (sum(keep) < 5L) { values(out) <- NA_real_; return(out) }
  nb_sub <- spdep::subset.nb(nb, keep)
  lw <- spdep::nb2listw(nb_sub, style = config$style, zero.policy = TRUE)
  z <- as.numeric(spdep::localG(x[valid][keep], listw = lw, zero.policy = TRUE))
  vals <- rep(NA_real_, length(valid))
  land <- rep(NA_real_, sum(valid))
  land[keep] <- z; vals[valid] <- land
  values(out) <- vals; out
}

gi_country <- function(r_layer, precomp) {
  x <- values(r_layer, mat = FALSE)
  x_land <- x[precomp$valid]
  z_land <- rep(NA_real_, sum(precomp$valid))
  for (s in precomp$country_structures) {
    xc <- x_land[s$mask_land]; keep <- is.finite(xc)
    if (sum(keep) < 5) next
    nb_sub <- spdep::subset.nb(s$nb, keep)
    lw <- spdep::nb2listw(nb_sub, style = config$style, zero.policy = TRUE)
    z <- as.numeric(spdep::localG(xc[keep], listw = lw, zero.policy = TRUE))
    z_land[s$indices_land[keep]] <- z
  }
  out <- rast(r_layer)
  vals <- rep(NA_real_, length(x))
  vals[precomp$valid] <- z_land
  values(out) <- vals; out
}

# ============================================================
# CLASSIFY + FDR (BH only)
# ============================================================

# Uncorrected: band by |z| against the z-cutoffs of gi_alphas
classify_z <- function(z) {
  cat <- rep(NA_real_, length(z)); f <- is.finite(z); cat[f] <- 0
  for (k in seq_along(Z_CUTS)) {
    cat[f & z >  Z_CUTS[k]] <-  k
    cat[f & z < -Z_CUTS[k]] <- -k
  }
  cat
}

# Corrected: band by the adjusted p-value (q) against gi_alphas directly.

classify_fdr <- function(q, z) {
  cat <- rep(NA_real_, length(q))
  f <- is.finite(q) & is.finite(z)
  cat[f] <- 0
  for (k in seq_along(config$gi_alphas)) {
    sel <- f & q <= config$gi_alphas[k] & z > 0
    cat[sel] <-  k
    sel <- f & q <= config$gi_alphas[k] & z < 0
    cat[sel] <- -k
  }
  cat
}

do_fdr <- function(z_land, p_land, mode, precomp = NULL) {
  res <- list()
  for (m in config$fdr_methods) {
    q <- rep(NA_real_, length(p_land))
    sig <- rep(NA_real_, length(p_land))
    cat_fdr <- rep(NA_real_, length(z_land))

    if (mode == "ssa") {
      fin <- is.finite(p_land)
      q[fin] <- p.adjust(p_land[fin], method = m)
      sig[fin] <- as.numeric(q[fin] <= config$fdr_alpha)
      cat_fdr <- classify_fdr(q, z_land)
    } else {
      for (s in precomp$country_structures) {
        idx <- which(precomp$country_vals_land == s$country_id & is.finite(p_land))
        if (!length(idx)) next
        q_c <- p.adjust(p_land[idx], method = m)
        q[idx] <- q_c
        sig[idx] <- as.numeric(q_c <= config$fdr_alpha)
        cat_fdr[idx] <- classify_fdr(q_c, z_land[idx])
      }
    }
    res[[m]] <- list(q = q, sig = sig, cat = cat_fdr)
  }
  res
}

# ============================================================
# FULL Gi* FOR ONE LAYER
# ============================================================

gi_full <- function(r_layer, nb_or_precomp, valid_v, mode) {
  if (mode == "ssa") {
    z_rast <- gi_ssa(r_layer, nb_or_precomp, valid_v)
    z_land <- values(z_rast, mat = FALSE)[valid_v]
    vm <- valid_v
  } else {
    z_rast <- gi_country(r_layer, valid_v)
    z_land <- values(z_rast, mat = FALSE)[valid_v$valid]
    vm <- valid_v$valid
  }
  fin <- is.finite(z_land)
  p_land <- rep(NA_real_, length(z_land))
  p_land[fin] <- 2 * pnorm(-abs(z_land[fin]))
  fdr <- do_fdr(z_land, p_land, mode, if (mode == "country") valid_v else NULL)
  cat_land <- classify_z(z_land)

  to_rast <- function(v, nm) {
    r <- rast(z_rast); full <- rep(NA_real_, ncell(r))
    full[vm] <- v; values(r) <- full; names(r) <- nm; r
  }

  out <- list(z = z_rast, p = to_rast(p_land, "p"), cat = to_rast(cat_land, "cat"))
  for (m in config$fdr_methods) {
    out[[paste0("q_", m)]]       <- to_rast(fdr[[m]]$q,   paste0("q_", m))
    out[[paste0("sig_fdr_", m)]] <- to_rast(fdr[[m]]$sig, paste0("sig_", m))
    out[[paste0("cat_fdr_", m)]] <- to_rast(fdr[[m]]$cat, paste0("cat_", m))
  }
  out
}

# ============================================================
# Process all layers
# ============================================================

gi_all_layers <- function(r_stack, layer_names, nb_or_precomp, valid_v, mode) {
  present <- intersect(layer_names, names(r_stack))
  missing <- setdiff(layer_names, names(r_stack))
  if (!length(present)) stop("The input stack holds no indicator or prevalence layers; rerun step 16.")
  if (length(missing))
    message(sprintf("  Warning: %d missing: %s", length(missing), paste(missing, collapse = ", ")))
  message(sprintf("  Processing %d layers (%s)...", length(present), toupper(mode)))

  keys <- c("z", "p", "cat")
  for (m in config$fdr_methods)
    keys <- c(keys, paste0("q_", m), paste0("sig_fdr_", m), paste0("cat_fdr_", m))
  coll <- setNames(lapply(keys, function(k) vector("list", length(present))), keys)

  for (i in seq_along(present)) {
    nm <- present[i]
    res <- gi_full(r_stack[[nm]], nb_or_precomp, valid_v, mode)
    names(res$z) <- paste0("GiZ_", nm)
    names(res$p) <- paste0("Gi_p_", nm)
    names(res$cat) <- paste0("Cat_", nm)
    coll$z[[i]] <- res$z; coll$p[[i]] <- res$p; coll$cat[[i]] <- res$cat
    for (m in config$fdr_methods) {
      names(res[[paste0("q_", m)]]) <- paste0("Gi_qFDR_", m, "_", nm)
      names(res[[paste0("sig_fdr_", m)]]) <- paste0("Gi_sigFDR_", m, "_", nm)
      names(res[[paste0("cat_fdr_", m)]]) <- paste0("CatFDR_", m, "_", nm)
      coll[[paste0("q_", m)]][[i]] <- res[[paste0("q_", m)]]
      coll[[paste0("sig_fdr_", m)]][[i]] <- res[[paste0("sig_fdr_", m)]]
      coll[[paste0("cat_fdr_", m)]][[i]] <- res[[paste0("cat_fdr_", m)]]
    }
    if (i %% 5 == 0) gc(verbose = FALSE)
  }
  lapply(coll, function(lst) do.call(c, lst))
}

# ============================================================
# Main processing
# ============================================================

message("  Locating the input file")

input_files <- list.files(config$in_dir, pattern = config$in_pattern, full.names = TRUE)
if (!length(input_files)) stop("No Gi* input stack found, run step 16 first. Looked in: ", config$in_dir)
all_file <- grep("ALL", input_files, value = TRUE)
if (!length(all_file)) stop("No combined-sex Gi* input stack; step 16 writes it into results/gistar/inputs/ready.")
in_file <- all_file[1]
message(sprintf("  Input: %s", basename(in_file)))
base <- sub(config$in_pattern, "", basename(in_file))

t_start <- Sys.time()

# --- Load ---
r <- rast(in_file)
message(sprintf("  %d layers", nlyr(r)))

# --- PREV summaries ---
draw_idx <- grep(config$prev_draw_pattern, names(r))
message(sprintf("  %d draws found", length(draw_idx)))
if (!length(draw_idx)) stop("The input stack holds no prevalence draw layers; rerun step 16.")

prev_draws <- r[[draw_idx]]
t0 <- Sys.time()
prev_mean  <- app(prev_draws, fun = "mean", na.rm = TRUE); names(prev_mean) <- "PREVpct_mean"
prev_lower <- app(prev_draws, fun = function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE))
names(prev_lower) <- "PREVpct_lower"
prev_upper <- app(prev_draws, fun = function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE))
names(prev_upper) <- "PREVpct_upper"
message(sprintf("  Done in %.1f sec", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
rm(prev_draws); gc(FALSE)

# --- Build analysis stack ---
ewe_present <- intersect(config$ewe_layers, names(r))
analysis_stack <- c(r[[ewe_present]], prev_mean, prev_lower, prev_upper, r$country)
all_layers <- c(ewe_present, "PREVpct_mean", "PREVpct_lower", "PREVpct_upper")
message(sprintf("  %d data layers + country", length(all_layers)))

country_rast <- r$country

# Country code -> name, read from the layer's own level table.
country_levels <- get_country_levels(country_rast)
message(sprintf("  Country lookup: %d countries from raster factor levels",
                nrow(country_levels)))

id_to_name <- function(cid) country_name(cid, country_levels)

rm(r, prev_mean, prev_lower, prev_upper); gc(FALSE)

# === SSA-WIDE ===
nb_ssa <- build_nb_ssa(analysis_stack)
res_ssa <- gi_all_layers(analysis_stack, all_layers, nb_ssa$nb, nb_ssa$valid, "ssa")

main <- c(res_ssa$z, res_ssa$p)
for (m in config$fdr_methods)
  main <- c(main, res_ssa[[paste0("q_", m)]], res_ssa[[paste0("sig_fdr_", m)]])
out_ssa_res <- file.path(config$out_dir, paste0(base, "_SSA", config$out_suffix_results))
writeRaster(main, out_ssa_res, overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"))
message(sprintf("  Written: %s (%d layers)", basename(out_ssa_res), nlyr(main)))

cats <- res_ssa$cat
for (m in config$fdr_methods) cats <- c(cats, res_ssa[[paste0("cat_fdr_", m)]])
out_ssa_cat <- file.path(config$out_dir, paste0(base, "_SSA", config$out_suffix_categorical))
writeRaster(cats, out_ssa_cat, overwrite = TRUE, datatype = "INT2S",
            gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "NUM_THREADS=ALL_CPUS"))
message(sprintf("  Written: %s (%d layers)", basename(out_ssa_cat), nlyr(cats)))
rm(res_ssa, main, cats, nb_ssa); gc(FALSE)

# === COUNTRY ===
precomp <- build_nb_country(analysis_stack)
res_ctry <- gi_all_layers(analysis_stack, all_layers, NULL, precomp, "country")

main <- c(res_ctry$z, res_ctry$p)
for (m in config$fdr_methods)
  main <- c(main, res_ctry[[paste0("q_", m)]], res_ctry[[paste0("sig_fdr_", m)]])
out_ctry_res <- file.path(config$out_dir, paste0(base, "_COUNTRY", config$out_suffix_results))
writeRaster(main, out_ctry_res, overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"))
message(sprintf("  Written: %s (%d layers)", basename(out_ctry_res), nlyr(main)))

cats <- res_ctry$cat
for (m in config$fdr_methods) cats <- c(cats, res_ctry[[paste0("cat_fdr_", m)]])
out_ctry_cat <- file.path(config$out_dir, paste0(base, "_COUNTRY", config$out_suffix_categorical))
writeRaster(cats, out_ctry_cat, overwrite = TRUE, datatype = "INT2S",
            gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "NUM_THREADS=ALL_CPUS"))
message(sprintf("  Written: %s (%d layers)", basename(out_ctry_cat), nlyr(cats)))
rm(res_ctry, main, cats, precomp, analysis_stack); gc(FALSE)

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
message(sprintf("\n  Analysis complete in %.1f minutes", elapsed))

# ============================================================
# PLHIV in hotspots and coldspots
# ============================================================

plhiv_med <- NULL
plhiv_mat <- NULL  # ncells x 1000 matrix

if (!file.exists(config$plhiv_ui_path)) {
  message("  Warning: PLHIV interval file not found, skipping population counts")
  plhiv_med <- NULL
} else {
  ui_r <- rast(config$plhiv_ui_path)
  plhiv_idx <- grep("PLHIV_mean|PLHIV_median|PLHIV_med", names(ui_r))
  if (length(plhiv_idx) > 0) plhiv_med <- values(ui_r[[plhiv_idx[1]]], mat = FALSE)
  rm(ui_r)
}

if (file.exists(config$plhiv_draws_path)) {
  message("  Loading 1000 PLHIV draws for CI...")
  draws_r <- rast(config$plhiv_draws_path)
  draw_idx <- grep("^PLHIV_d[0-9]", names(draws_r))
  if (length(draw_idx) > 0) {
    plhiv_mat <- values(draws_r[[draw_idx]], mat = TRUE)
    plhiv_mat[!is.finite(plhiv_mat)] <- 0
    message(sprintf("  PLHIV draws: %d cells x %d draws", nrow(plhiv_mat), ncol(plhiv_mat)))
    # Use draw mean as plhiv_med if not already loaded
    if (is.null(plhiv_med)) plhiv_med <- rowMeans(plhiv_mat, na.rm = TRUE)
  }
  rm(draws_r); gc(FALSE)
} else {
  message("  Warning: PLHIV draws file not found, reporting the mean without intervals")
}

country_v <- values(country_rast, mat = FALSE)
if (!is.null(plhiv_med)) {
  message(sprintf("  PLHIV loaded: %s total, %s cells",
                  format(sum(plhiv_med, na.rm = TRUE), big.mark = ","),
                  format(sum(!is.na(plhiv_med)), big.mark = ",")))
}

fmt_ci <- function(med, lo, hi) {
  f <- function(x) {
    if (is.na(x)) return("NA")
    if (abs(x) >= 1e6) return(sprintf("%.1fM", x/1e6))
    if (abs(x) >= 1e3) return(sprintf("%.0fK", x/1e3))
    sprintf("%.0f", x)
  }
  if (is.na(lo)) return(f(med))
  paste0(f(med), " [", f(lo), "\u2013", f(hi), "]")
}

if (!is.null(plhiv_med)) {

  plhiv_rows <- list()

  for (ref in c("SSA", "COUNTRY")) {
    cat_fp <- file.path(config$out_dir,
                        paste0(base, "_", ref, config$out_suffix_categorical))
    if (!file.exists(cat_fp)) next
    r_cat <- rast(cat_fp)

    # Get BH categorical layers
    bh_lyrs <- grep("^CatFDR_BH_", names(r_cat), value = TRUE)

    for (cl in bh_lyrs) {
      var <- sub("^CatFDR_BH_", "", cl)
      cat_v <- values(r_cat[[cl]], mat = FALSE)

      valid <- !is.na(cat_v) & !is.na(plhiv_med)
      # "Significant" == survives BH at config$fdr_alpha (code magnitude SIG_CODE)
      hot_idx  <- valid & cat_v >=  SIG_CODE
      cold_idx <- valid & cat_v <= -SIG_CODE
      ns_idx   <- valid & abs(cat_v) < SIG_CODE

      # SSA with CI
      ci_total <- count_plhiv_ci(valid)
      ci_hot   <- count_plhiv_ci(hot_idx)
      ci_cold  <- count_plhiv_ci(cold_idx)
      pci_h <- pct_plhiv_ci(hot_idx, valid)
      pci_c <- pct_plhiv_ci(cold_idx, valid)

      plhiv_rows[[length(plhiv_rows) + 1]] <- data.frame(
        mode = ref, variable = var, label = get_label(var),
        region = "SSA", correction = "BH", alpha = config$fdr_alpha,
        n_total = sum(valid), n_hotspot = sum(hot_idx), n_coldspot = sum(cold_idx),
        plhiv_total_mean = ci_total$plhiv_mean, plhiv_total_lower = ci_total$plhiv_lower, plhiv_total_upper = ci_total$plhiv_upper,
        plhiv_hotspot_mean = ci_hot$plhiv_mean, plhiv_hotspot_lower = ci_hot$plhiv_lower, plhiv_hotspot_upper = ci_hot$plhiv_upper,
        plhiv_coldspot_mean = ci_cold$plhiv_mean, plhiv_coldspot_lower = ci_cold$plhiv_lower, plhiv_coldspot_upper = ci_cold$plhiv_upper,
        pct_hot_mean = pci_h$pct_mean, pct_hot_lower = pci_h$pct_lower, pct_hot_upper = pci_h$pct_upper,
        pct_cold_mean = pci_c$pct_mean, pct_cold_lower = pci_c$pct_lower, pct_cold_upper = pci_c$pct_upper,
        stringsAsFactors = FALSE)

      # Per-country
      country_ids <- sort(unique(country_v[valid & !is.na(country_v)]))
      for (cid in country_ids) {
        c_mask <- valid & country_v == cid
        c_hot  <- c_mask & cat_v >=  SIG_CODE
        c_cold <- c_mask & cat_v <= -SIG_CODE

        ci_ct <- count_plhiv_ci(c_mask)
        if (ci_ct$plhiv_mean < 1) next
        ci_ch <- count_plhiv_ci(c_hot)
        ci_cc <- count_plhiv_ci(c_cold)
        pci_ch2 <- pct_plhiv_ci(c_hot, c_mask)
        pci_cc2 <- pct_plhiv_ci(c_cold, c_mask)

        plhiv_rows[[length(plhiv_rows) + 1]] <- data.frame(
          mode = ref, variable = var, label = get_label(var),
          region = id_to_name(cid), correction = "BH", alpha = config$fdr_alpha,
          n_total = sum(c_mask), n_hotspot = sum(c_hot), n_coldspot = sum(c_cold),
          plhiv_total_mean = ci_ct$plhiv_mean, plhiv_total_lower = ci_ct$plhiv_lower, plhiv_total_upper = ci_ct$plhiv_upper,
          plhiv_hotspot_mean = ci_ch$plhiv_mean, plhiv_hotspot_lower = ci_ch$plhiv_lower, plhiv_hotspot_upper = ci_ch$plhiv_upper,
          plhiv_coldspot_mean = ci_cc$plhiv_mean, plhiv_coldspot_lower = ci_cc$plhiv_lower, plhiv_coldspot_upper = ci_cc$plhiv_upper,
          pct_hot_mean = pci_ch2$pct_mean, pct_hot_lower = pci_ch2$pct_lower, pct_hot_upper = pci_ch2$pct_upper,
          pct_cold_mean = pci_cc2$pct_mean, pct_cold_lower = pci_cc2$pct_lower, pct_cold_upper = pci_cc2$pct_upper,
          stringsAsFactors = FALSE)
      }
    }
    rm(r_cat); gc(FALSE)
  }

  plhiv_df <- do.call(rbind, plhiv_rows)
  plhiv_csv <- file.path(config$out_dir, "PLHIV_in_hotspots_coldspots.csv")
  write.csv(plhiv_df, plhiv_csv, row.names = FALSE)
  message(sprintf("\n  PLHIV CSV saved: %s (%d rows)", basename(plhiv_csv), nrow(plhiv_df)))

  ssa_sub <- plhiv_df[plhiv_df$region == "SSA", ]
  message(sprintf("\n  PLHIV IN HOTSPOTS - SSA (BH corrected, q <= %.2f, mean [95%% UI]):",
                  config$fdr_alpha))
  for (i in seq_len(nrow(ssa_sub))) {
    message(sprintf("  %-30s  Hotspot: %s  (%.1f%%)",
                    ssa_sub$label[i],
                    fmt_ci(ssa_sub$plhiv_hotspot_mean[i], ssa_sub$plhiv_hotspot_lower[i], ssa_sub$plhiv_hotspot_upper[i]),
                    ssa_sub$pct_hot_mean[i]))
  }
}

# ============================================================
# Maps
# ============================================================

message("  Maps")

viz <- list(
  out_dir  = config$out_dir,
  plot_dir = DIR_GI_PLOTS,
  width    = 12,    # 2-panel (wider per panel)
  height   = 7,
  z_width  = 10,
  z_height = 8
)

dir.create(viz$plot_dir, showWarnings = FALSE, recursive = TRUE)
for (d in c("SSA", "COUNTRY", "prevalence_estimates"))
  dir.create(file.path(viz$plot_dir, d), showWarnings = FALSE, recursive = TRUE)

# --- Load GADM boundaries ---
borders_v <- NULL
if (!is.null(config$gadm_path) && file.exists(config$gadm_path)) {
  borders_sf <- sf::st_read(config$gadm_path, quiet = TRUE)
  borders_v <- vect(borders_sf)
  message(sprintf("  GADM boundaries loaded: %d country polygons", nrow(borders_v)))
}

# --- Colour palettes ---
col_div <- function(n = 100) colorRampPalette(c("blue", "lightblue", "white", "pink", "red"))(n)

# Hotspot colours: codes -3..+3
col_hot_full <- c(
  "#08519c",   # -3  (cold 99%)
  "#3182bd",   # -2  (cold 95%)
  "#9ecae1",   # -1  (cold 90%)
  "#F5F0E1",   # 0   (not significant)
  "#fcae91",   # +1  (hot 90%)
  "#fb6a4a",   # +2  (hot 95%)
  "#cb181d"    # +3  (hot 99%)
)
hot_breaks <- c(-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5)

# Uncorrected panels are thresholded on raw p; BH panels on the adjusted q.
hot_legend <- c("Cold (p<0.01)", "Cold (p<0.05)", "Cold (p<0.10)", "Not significant",
                "Hot (p<0.10)", "Hot (p<0.05)", "Hot (p<0.01)")
hot_legend_bh <- c("Cold (q<0.01)", "Cold (q<0.05)", "Cold (q<0.10)", "Not significant",
                   "Hot (q<0.10)", "Hot (q<0.05)", "Hot (q<0.01)")
hot_legend_gen <- c("Cold (99%)", "Cold (95%)", "Cold (90%)", "Not significant",
                    "Hot (90%)", "Hot (95%)", "Hot (99%)")

# Hotspot-only: reclassify coldspots as not significant
col_hot_only <- c(
  "#F5F0E1",   # -3  coldspots are drawn as not significant on this variant
  "#F5F0E1",   # -2
  "#F5F0E1",   # -1
  "#F5F0E1",   # 0   (not significant)
  "#fcae91",   # +1  (hot 90%)
  "#fb6a4a",   # +2  (hot 95%)
  "#cb181d"    # +3  (hot 99%)
)
hot_only_legend    <- c("Not a hotspot", "A hotspot (p<0.10)", "A hotspot (p<0.05)", "A hotspot (p<0.01)")
hot_only_legend_bh <- c("Not a hotspot", "A hotspot (q<0.10)", "A hotspot (q<0.05)", "A hotspot (q<0.01)")
hot_only_cols   <- c("#F5F0E1", "#fcae91", "#fb6a4a", "#cb181d")

# --- Draw boundary helper ---
add_boundary <- function(ref) {
  if (!is.null(borders_v)) {
    lwd <- if (ref == "COUNTRY") 0.8 else 0.5
    lines(borders_v, col = "black", lwd = lwd)
  }
}

# --- Find output files ---
res_esc <- gsub("\\.", "\\\\.", config$out_suffix_results)
cat_esc <- gsub("\\.", "\\\\.", config$out_suffix_categorical)

cat_files <- list.files(viz$out_dir, pattern = paste0(cat_esc, "$"), full.names = TRUE)
res_files <- list.files(viz$out_dir, pattern = paste0(res_esc, "$"), full.names = TRUE)

cat_ssa <- grep("_SSA", cat_files, value = TRUE)
cat_cty <- grep("_COUNTRY", cat_files, value = TRUE)
res_ssa <- grep("_SSA", res_files, value = TRUE)
res_cty <- grep("_COUNTRY", res_files, value = TRUE)

message(sprintf("  Found: %d categorical, %d results files",
                length(cat_files), length(res_files)))

# ============================================================
# 1. SIGNIFICANCE COUNTS TABLE
# ============================================================


sig_counts_all <- list()

for (fp in c(cat_ssa, cat_cty)) {
  ref <- if (grepl("_SSA", basename(fp))) "SSA" else "COUNTRY"
  r <- rast(fp)

  cat_lyrs <- grep("^Cat_[^F]", names(r), value = TRUE)

  for (cl in cat_lyrs) {
    var <- sub("^Cat_", "", cl)
    v <- values(r[[cl]], mat = FALSE)
    v <- v[!is.na(v)]
    n_total <- length(v)

    n_hot_90 <- sum(v >= 1); n_cold_90 <- sum(v <= -1)
    n_hot_95 <- sum(v >= 2); n_cold_95 <- sum(v <= -2)
    n_hot_99 <- sum(v >= 3); n_cold_99 <- sum(v <= -3)

    row_base <- data.frame(mode = ref, variable = var, n_land = n_total,
                           stringsAsFactors = FALSE)

    sig_counts_all[[length(sig_counts_all) + 1]] <- cbind(row_base, data.frame(
      correction = "Uncorrected",
      n_hot_90 = n_hot_90, n_hot_95 = n_hot_95, n_hot_99 = n_hot_99,
      n_cold_90 = n_cold_90, n_cold_95 = n_cold_95, n_cold_99 = n_cold_99,
      pct_sig_95 = round(100 * (n_hot_95 + n_cold_95) / n_total, 2)
    ))

    # BH
    fdr_lyr <- paste0("CatFDR_BH_", var)
    if (fdr_lyr %in% names(r)) {
      vf <- values(r[[fdr_lyr]], mat = FALSE)
      vf <- vf[!is.na(vf)]

      sig_counts_all[[length(sig_counts_all) + 1]] <- cbind(row_base, data.frame(
        correction = "BH",
        n_hot_90 = sum(vf >= 1), n_hot_95 = sum(vf >= 2), n_hot_99 = sum(vf >= 3),
        n_cold_90 = sum(vf <= -1), n_cold_95 = sum(vf <= -2), n_cold_99 = sum(vf <= -3),
        pct_sig_95 = round(100 * (sum(abs(vf) >= 2)) / n_total, 2)
      ))
    }
  }
  rm(r); gc(FALSE)
}

sig_counts_df <- do.call(rbind, sig_counts_all)
sig_csv <- file.path(viz$plot_dir, "significance_counts.csv")
write.csv(sig_counts_df, sig_csv, row.names = FALSE)
message(sprintf("  Saved: %s", basename(sig_csv)))

# Print summary
message("\n  Significance before and after correction, 95% level (p < 0.05 against q < 0.05):")
message(sprintf("  %-12s %-30s %8s %8s %6s",
                "Mode", "Variable", "Uncorr", "BH", "Retain%"))
for (var in unique(sig_counts_df$variable)) {
  for (ref in c("SSA", "COUNTRY")) {
    sub <- sig_counts_df[sig_counts_df$variable == var & sig_counts_df$mode == ref, ]
    if (nrow(sub) == 0) next
    unc <- sub$n_hot_95[sub$correction == "Uncorrected"] +
      sub$n_cold_95[sub$correction == "Uncorrected"]
    hlm <- sub$n_hot_95[sub$correction == "BH"] +
      sub$n_cold_95[sub$correction == "BH"]
    pct <- if (length(hlm) && length(unc) && unc > 0) round(100 * hlm / unc, 1) else NA
    message(sprintf("  %-12s %-30s %8s %8s %5.1f%%",
                    ref, get_label(var),
                    format(unc, big.mark = ","),
                    format(hlm, big.mark = ","),
                    ifelse(is.na(pct), NA, pct)))
  }
}

# ============================================================
# 2. Z-SCORE MAP + 2-PANEL CATEGORICAL (Uncorrected | BH)
# ============================================================


plot_var <- function(res_fp, cat_fp, ref) {
  id <- sub(paste0(cat_esc, "$"), "", basename(cat_fp))
  r_res <- rast(res_fp)
  r_cat <- rast(cat_fp)

  z_lyrs <- grep("^GiZ_", names(r_res), value = TRUE)

  for (zl in z_lyrs) {
    var <- sub("^GiZ_", "", zl)

    # --- Z-score map ---
    png(file.path(viz$plot_dir, ref, paste0(id, "_GiZ_", var, ".png")),
        width = viz$z_width, height = viz$z_height, units = "in", res = 600)
    par(mar = c(2, 2, 3, 4))
    plot(r_res[[zl]], col = col_div(100),
         main = "",
         axes = FALSE, zlim = c(-5, 5), cex.main = 1.3)
    add_boundary(ref)
    dev.off()

    # --- 2-panel categorical (Uncorrected | BH) ---
    cat_uncorr <- paste0("Cat_", var)
    cat_bh   <- paste0("CatFDR_BH_", var)

    panels <- c(cat_uncorr, cat_bh)
    labels <- c("Uncorrected", "BH corrected")
    panels_present <- panels %in% names(r_cat)

    if (sum(panels_present) >= 1) {
      panels <- panels[panels_present]
      labels <- labels[panels_present]

      png(file.path(viz$plot_dir, ref, paste0(id, "_Hotspot_", var, ".png")),
          width = 6 * length(panels), height = 7, units = "in", res = 600)
      np <- length(panels)
      layout(rbind(seq_len(np), rep(np + 1, np)), heights = c(1, 0.12))
      par(mar = c(2, 1, 3, 1), oma = c(0, 0, 3, 0))

      for (j in seq_along(panels)) {
        plot(r_cat[[panels[j]]], col = col_hot_full,
             main = labels[j], axes = FALSE,
             breaks = hot_breaks, legend = FALSE,
             colNA = "white", cex.main = 1.3)
        add_boundary(ref)
      }

      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = hot_legend_gen, fill = col_hot_full,
             ncol = 4, cex = 1.1, bty = "n", xpd = TRUE,
             text.width = NA, x.intersp = 0.8)

      dev.off()
    }

    # --- Standalone BH-only map ---
    if (cat_bh %in% names(r_cat)) {
      jpeg(file.path(viz$plot_dir, ref, paste0(id, "_BH_only_", var, ".jpeg")),
           width = 10, height = 9, units = "in", res = 600, quality = 95)
      layout(matrix(1:2, 2, 1), heights = c(1, 0.08))
      par(mar = c(1, 2, 3, 2))
      plot(r_cat[[cat_bh]], col = col_hot_full,
           main = "",
           axes = FALSE, breaks = hot_breaks, legend = FALSE,
           colNA = "white", cex.main = 1.2)
      add_boundary(ref)
      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = hot_legend_bh, fill = col_hot_full,
             cex = 0.75, bty = "n", ncol = 4)
      dev.off()

      # --- Hotspot-only version ---
      jpeg(file.path(viz$plot_dir, ref, paste0(id, "_BH_hotonly_", var, ".jpeg")),
           width = 10, height = 9, units = "in", res = 600, quality = 95)
      layout(matrix(1:2, 2, 1), heights = c(1, 0.08))
      par(mar = c(1, 2, 3, 2))
      plot(r_cat[[cat_bh]], col = col_hot_only,
           main = "",
           axes = FALSE, breaks = hot_breaks, legend = FALSE,
           colNA = "white", cex.main = 1.2)
      add_boundary(ref)
      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = hot_only_legend_bh, fill = hot_only_cols,
             cex = 0.85, bty = "n", ncol = 4)
      dev.off()
    }
  }

  rm(r_res, r_cat); gc(FALSE)
}

# -- BH Panel Maps: 3 main + 6 sensitivity --
make_bh_panel <- function(cat_fp, ref, vars, panel_tag, panel_title) {
  r_cat <- rast(cat_fp)
  nv <- length(vars)
  present <- character()
  for (v in vars) {
    nm <- paste0("CatFDR_BH_", v)
    if (nm %in% names(r_cat)) present <- c(present, v)
  }
  if (length(present) < 2) { rm(r_cat); return() }
  nr <- ceiling(length(present) / 3); nc <- min(3, length(present))
  jpeg(file.path(viz$plot_dir, ref, paste0("BH_panel_", panel_tag, "_", ref, ".jpeg")),
       width = 5 * nc, height = 4 * nr + 3, units = "in", res = 600, quality = 95)
  layout(rbind(matrix(seq_len(nr * nc), nr, nc, byrow = TRUE), rep(nr * nc + 1, nc)), heights = c(rep(1, nr), 0.12))
  par(mar = c(1, 1, 2.5, 1), oma = c(0, 0, 3, 0))
  abc <- letters[seq_along(present)]
  for (i in seq_along(present)) {
    nm <- paste0("CatFDR_BH_", present[i])
    plot(r_cat[[nm]], col = col_hot_full, breaks = hot_breaks,
         main = paste0(abc[i], ") ", get_label(present[i])),
         axes = FALSE, legend = FALSE, colNA = "white", cex.main = 1.0)
    add_boundary(ref)
  }
  remainder <- (nc * nr) - length(present)
  if (remainder > 0) for (r in seq_len(remainder)) plot.new()
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = hot_legend_bh, fill = col_hot_full,
         ncol = 4, cex = 0.85, bty = "n", x.intersp = 0.5)
  dev.off()

  # --- Hotspot-only panel ---
  jpeg(file.path(viz$plot_dir, ref, paste0("BH_panel_hotonly_", panel_tag, "_", ref, ".jpeg")),
       width = 5 * nc, height = 4 * nr + 3, units = "in", res = 600, quality = 95)
  layout(rbind(matrix(seq_len(nr * nc), nr, nc, byrow = TRUE), rep(nr * nc + 1, nc)), heights = c(rep(1, nr), 0.12))
  par(mar = c(1, 1, 2.5, 1), oma = c(0, 0, 3, 0))
  for (i in seq_along(present)) {
    nm <- paste0("CatFDR_BH_", present[i])
    plot(r_cat[[nm]], col = col_hot_only, breaks = hot_breaks,
         main = paste0(abc[i], ") ", get_label(present[i])),
         axes = FALSE, legend = FALSE, colNA = "white", cex.main = 1.0)
    add_boundary(ref)
  }
  if (remainder > 0) for (r in seq_len(remainder)) plot.new()
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = hot_only_legend_bh, fill = hot_only_cols,
         ncol = 4, cex = 0.85, bty = "n", x.intersp = 0.5)
  dev.off()
  rm(r_cat); gc(FALSE)
  message(sprintf("  BH panel %s [%s] saved", panel_tag, ref))
}

main_3 <- config$main_vars
sens_6 <- setdiff(config$ewe_layers, config$main_vars)
sens_6 <- sens_6[!grepl("PREV", sens_6)]

for (ref_mode in c("SSA", "COUNTRY")) {
  cat_fps <- if (ref_mode == "SSA") cat_ssa else cat_cty
  if (length(cat_fps) == 0) next
  make_bh_panel(cat_fps[1], ref_mode, main_3, "main3", "Main EWE Indicators")
  make_bh_panel(cat_fps[1], ref_mode, sens_6, "sens6", "Sensitivity EWE Indicators")
}

# ------------------------------------------------------------
# PLHIV in hotspots and coldspots at three significance levels
# ------------------------------------------------------------
# Thresholding is done here on the stored q and p values rather than on the
# category rasters, so the three levels come from one source and a different
# level could be added without recomputing anything.

message("  PLHIV in Gi* hotspots and coldspots, by threshold")

# Significance levels (two-sided alpha), NOT category codes
gi_thresholds <- c("90%" = 0.10, "95%" = 0.05, "99%" = 0.01)

# Which p-value columns to threshold: BH-corrected q, and raw uncorrected p
gi_corrections <- c("BH", "uncorrected")

if (!is.null(plhiv_mat)) {

  gi_plhiv_rows <- list()

  for (ref_mode in c("SSA", "COUNTRY")) {
    res_fps <- if (ref_mode == "SSA") res_ssa else res_cty
    if (length(res_fps) == 0) next

    # q, p and z are in the results raster, not the categorical raster
    r_res <- rast(res_fps[1])
    z_lyrs <- grep("^GiZ_", names(r_res), value = TRUE)

    for (zl in z_lyrs) {
      var <- sub("^GiZ_", "", zl)

      z_v <- values(r_res[[zl]], mat = FALSE)

      for (corr in gi_corrections) {
        pl <- if (corr == "BH") paste0("Gi_qFDR_BH_", var) else paste0("Gi_p_", var)
        if (!pl %in% names(r_res)) next
        p_v <- values(r_res[[pl]], mat = FALSE)

        valid <- is.finite(z_v) & is.finite(p_v) & !is.na(plhiv_med)
        if (!any(valid)) next

        # Denominators do not depend on alpha - compute once per (var, correction)
        ci_total <- count_plhiv_ci(valid)
        cids <- sort(unique(country_v[valid & !is.na(country_v) & country_v > 0]))
        c_valid_list <- setNames(lapply(cids, function(cid)
          valid & !is.na(country_v) & country_v == cid), as.character(cids))
        ci_ct_list <- lapply(c_valid_list, count_plhiv_ci)

        for (thr_name in names(gi_thresholds)) {
          alpha <- gi_thresholds[[thr_name]]

          hot_mask  <- valid & p_v <= alpha & z_v > 0
          cold_mask <- valid & p_v <= alpha & z_v < 0

          ci_hot   <- count_plhiv_ci(hot_mask)
          ci_cold  <- count_plhiv_ci(cold_mask)

          pci_hot  <- pct_plhiv_ci(hot_mask, valid)
          pci_cold <- pct_plhiv_ci(cold_mask, valid)

          # SSA row
          gi_plhiv_rows[[length(gi_plhiv_rows) + 1]] <- data.frame(
            mode = ref_mode, variable = var, label = get_label(var),
            correction = corr, threshold = thr_name, alpha = alpha, region = "SSA",
            n_total = sum(valid), n_hotspot = sum(hot_mask), n_coldspot = sum(cold_mask),
            plhiv_total_mean = ci_total$plhiv_mean, plhiv_total_lower = ci_total$plhiv_lower, plhiv_total_upper = ci_total$plhiv_upper,
            plhiv_hot_mean = ci_hot$plhiv_mean, plhiv_hot_lower = ci_hot$plhiv_lower, plhiv_hot_upper = ci_hot$plhiv_upper,
            plhiv_cold_mean = ci_cold$plhiv_mean, plhiv_cold_lower = ci_cold$plhiv_lower, plhiv_cold_upper = ci_cold$plhiv_upper,
            pct_hot_mean = pci_hot$pct_mean, pct_hot_lower = pci_hot$pct_lower, pct_hot_upper = pci_hot$pct_upper,
            pct_cold_mean = pci_cold$pct_mean, pct_cold_lower = pci_cold$pct_lower, pct_cold_upper = pci_cold$pct_upper,
            stringsAsFactors = FALSE)

          # Per-country rows
          for (cid in cids) {
            c_valid <- c_valid_list[[as.character(cid)]]
            ci_ct   <- ci_ct_list[[as.character(cid)]]
            if (!any(c_valid)) next
            c_hot  <- c_valid & p_v <= alpha & z_v > 0
            c_cold <- c_valid & p_v <= alpha & z_v < 0

            ci_ch  <- count_plhiv_ci(c_hot)
            ci_cc  <- count_plhiv_ci(c_cold)
            pci_ch <- pct_plhiv_ci(c_hot, c_valid)
            pci_cc <- pct_plhiv_ci(c_cold, c_valid)

            gi_plhiv_rows[[length(gi_plhiv_rows) + 1]] <- data.frame(
              mode = ref_mode, variable = var, label = get_label(var),
              correction = corr, threshold = thr_name, alpha = alpha, region = id_to_name(cid),
              n_total = sum(c_valid), n_hotspot = sum(c_hot), n_coldspot = sum(c_cold),
              plhiv_total_mean = ci_ct$plhiv_mean, plhiv_total_lower = ci_ct$plhiv_lower, plhiv_total_upper = ci_ct$plhiv_upper,
              plhiv_hot_mean = ci_ch$plhiv_mean, plhiv_hot_lower = ci_ch$plhiv_lower, plhiv_hot_upper = ci_ch$plhiv_upper,
              plhiv_cold_mean = ci_cc$plhiv_mean, plhiv_cold_lower = ci_cc$plhiv_lower, plhiv_cold_upper = ci_cc$plhiv_upper,
              pct_hot_mean = pci_ch$pct_mean, pct_hot_lower = pci_ch$pct_lower, pct_hot_upper = pci_ch$pct_upper,
              pct_cold_mean = pci_cc$pct_mean, pct_cold_lower = pci_cc$pct_lower, pct_cold_upper = pci_cc$pct_upper,
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    rm(r_res); gc(FALSE)
  }

  gi_plhiv_df <- do.call(rbind, gi_plhiv_rows)
  gi_plhiv_csv <- file.path(config$out_dir, "PLHIV_Gi_multithreshold_hotcold.csv")
  write.csv(gi_plhiv_df, gi_plhiv_csv, row.names = FALSE)
  message(sprintf("  Saved: %s (%d rows)", basename(gi_plhiv_csv), nrow(gi_plhiv_df)))

  chk <- gi_plhiv_df[gi_plhiv_df$region == "SSA", ]
  for (k in split(chk, list(chk$mode, chk$variable, chk$correction), drop = TRUE)) {
    k <- k[order(-k$alpha), ]   # 0.10, 0.05, 0.01
    if (any(diff(k$n_hotspot) > 0) || any(diff(k$n_coldspot) > 0))
      message(sprintf("  Warning: non-nested thresholds for %s / %s / %s",
                      k$mode[1], k$variable[1], k$correction[1]))
    if (length(unique(k$n_hotspot)) == 1 && k$n_hotspot[1] > 0)
      message(sprintf("  NOTE: identical hotspot counts across 90/95/99 for %s / %s / %s",
                      k$mode[1], k$variable[1], k$correction[1]))
  }

  # Print SSA summary (BH only, to keep it readable)
  ssa_gi <- gi_plhiv_df[gi_plhiv_df$region == "SSA" & gi_plhiv_df$correction == "BH", ]
  message("\n  Gi* PLHIV - SSA, BH-corrected (mean [95% UI]):")
  for (i in seq_len(nrow(ssa_gi))) {
    message(sprintf("  %-25s %5s  Hot: %s  Cold: %s",
                    ssa_gi$label[i], ssa_gi$threshold[i],
                    fmt_ci(ssa_gi$plhiv_hot_mean[i], ssa_gi$plhiv_hot_lower[i], ssa_gi$plhiv_hot_upper[i]),
                    fmt_ci(ssa_gi$plhiv_cold_mean[i], ssa_gi$plhiv_cold_lower[i], ssa_gi$plhiv_cold_upper[i])))
  }
}

# Plot SSA
if (length(res_ssa) && length(cat_ssa)) {
  for (i in seq_along(res_ssa)) plot_var(res_ssa[i], cat_ssa[i], "SSA")
}

# Plot COUNTRY
if (length(res_cty) && length(cat_cty)) {
  for (i in seq_along(res_cty)) plot_var(res_cty[i], cat_cty[i], "COUNTRY")
}

# ============================================================
# 3. THE THREE HIV PREVALENCE ESTIMATES
# ------------------------------------------------------------
# The same Gi* analysis run on the mean prevalence surface and on the lower and
# upper bounds of its uncertainty interval so that the effect of
# prevalence uncertainty on the hotspot pattern can be read off directly.
# ============================================================


prev_vars <- c("PREVpct_mean", "PREVpct_lower", "PREVpct_upper")

for (ref in c("SSA", "COUNTRY")) {
  cat_fp <- if (ref == "SSA") cat_ssa[1] else cat_cty[1]
  if (is.null(cat_fp) || !file.exists(cat_fp)) next

  id <- sub(paste0(cat_esc, "$"), "", basename(cat_fp))
  r_cat <- rast(cat_fp)

  cat_names <- paste0("Cat_", prev_vars)
  if (!all(cat_names %in% names(r_cat))) {
    message(sprintf("    %s: Missing PREV categorical layers, skipping", ref))
    next
  }

  png(file.path(viz$plot_dir, "prevalence_estimates",
                paste0(id, "_PREV_mean_lower_upper.png")),
      width = 15, height = 7, units = "in", res = 600)
  layout(rbind(1:3, rep(4, 3)), heights = c(1, 0.12))
  par(mar = c(2, 1, 3, 1), oma = c(0, 0, 3, 0))

  prev_abc <- c("a) Mean", "b) Lower bound", "c) Upper bound")
  for (j in seq_along(prev_vars)) {
    cl <- paste0("Cat_", prev_vars[j])
    plot(r_cat[[cl]], col = col_hot_full, main = prev_abc[j],
         axes = FALSE, breaks = hot_breaks, legend = FALSE, colNA = "white",
         cex.main = 1.1)
    add_boundary(ref)
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = hot_legend, fill = col_hot_full,
         ncol = 4, cex = 1.1, bty = "n", xpd = TRUE)
  dev.off()

  # --- Hotspot-only version ---
  png(file.path(viz$plot_dir, "prevalence_estimates",
                paste0(id, "_PREV_mean_lower_upper_hotonly.png")),
      width = 15, height = 7, units = "in", res = 600)
  layout(rbind(1:3, rep(4, 3)), heights = c(1, 0.12))
  par(mar = c(2, 1, 3, 1), oma = c(0, 0, 3, 0))
  for (j in seq_along(prev_vars)) {
    cl <- paste0("Cat_", prev_vars[j])
    plot(r_cat[[cl]], col = col_hot_only, main = prev_abc[j],
         axes = FALSE, breaks = hot_breaks, legend = FALSE, colNA = "white",
         cex.main = 1.1)
    add_boundary(ref)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = hot_only_legend, fill = hot_only_cols,
         ncol = 4, cex = 1.1, bty = "n", xpd = TRUE)
  dev.off()

  # Quantify differences
  v_mean  <- values(r_cat[[cat_names[1]]], mat = FALSE)
  v_lower <- values(r_cat[[cat_names[2]]], mat = FALSE)
  v_upper <- values(r_cat[[cat_names[3]]], mat = FALSE)
  valid <- !is.na(v_mean) & !is.na(v_lower) & !is.na(v_upper)

  n_valid <- sum(valid)
  agree_all <- sum(v_mean[valid] == v_lower[valid] & v_mean[valid] == v_upper[valid])
  pct_agree <- round(100 * agree_all / n_valid, 2)
  robust_hot <- sum(v_mean[valid] >= 1 & v_lower[valid] >= 1 & v_upper[valid] >= 1)
  robust_cold <- sum(v_mean[valid] <= -1 & v_lower[valid] <= -1 & v_upper[valid] <= -1)

  message(sprintf("\n  PREV CI [%s]: %.1f%% agree | %s robust hot | %s robust cold",
                  ref, pct_agree,
                  format(robust_hot, big.mark = ","),
                  format(robust_cold, big.mark = ",")))

  rm(r_cat); gc(FALSE)
}


n_plots <- length(list.files(viz$plot_dir, pattern = "\\.png$", recursive = TRUE))
message(sprintf("  Plots: %d in %s", n_plots, viz$plot_dir))
message(sprintf("  Sig counts: %s", sig_csv))
if (exists("plhiv_csv")) message(sprintf("  PLHIV in hotspots: %s", plhiv_csv))

message("Done. Next: step 19.")
