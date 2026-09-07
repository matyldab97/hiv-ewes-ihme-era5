# ======================================================================
# Monte Carlo PLHIV draws
# ----------------------------------------------------------------------
# Uncertainty in the published IHME estimates is propagated with a three-
# component variance decomposition (Haeuser et al. 2022, BMC Medicine 20:488,
# doi 10.1186/s12916-022-02639-z, Supplementary Table S6):
#
#   Z1  spatial      one draw per spatial block, block size = the region's
#                    Matern range; shared across all 18 age-sex groups
#   Z2  age-sex      one 18-dimensional draw per region per iteration, with
#                    Kronecker correlation (sex x age); constant within region
#   Z3  country-age  two independent 9-dimensional draws per country per
#                    iteration, one per sex, AR(1) across age groups
#
#   Z_total[p, k] = sqrt(w1) Z1[p] + sqrt(w2) Z2[k] + sqrt(w3) Z3_c[k]
#
# with w1 + w2 + w3 = 1 the region-specific variance fractions. Z2 and Z3 are
# shared spatially and so do not cancel when fine cells are aggregated.
#
# PLHIV are drawn from a shifted lognormal, X = gamma + exp(mu + sigma Z),
# fitted to the published lower bound, mean and upper bound. Sampling is at
# the fine resolution, then aggregated to 0.25 degrees by area weighting.
# Population is held fixed across draws and prevalence is derived as
# 100 * PLHIV / population, rather than drawing the two independently. Drawing
# both would let a high PLHIV draw meet a low population draw in the same cell
# and produce prevalences that no combination of the published inputs supports.
# Point estimates are the mean of the draws; intervals are the 2.5th and
# 97.5th percentiles.
#
# Z1 uses discrete spatial blocks rather than a smooth Matern decay, Z3 is
# independent between sexes, and Table S6 medians are used in place of the
# full posterior.
#
# Input:  IHME gridded PLHIV and prevalence, GADM boundaries
# Output: results/montecarlo
# ======================================================================

source(file.path("R", "01_setup", "01_config.R"))
source(file.path("R", "01_setup", "03_lib_shared.R"))


library(terra)
library(sf)
library(exactextractr)
library(Matrix)
library(tools)

# 1000-draw stacks do not fit in memory, so every intermediate goes to disk.
terraOptions(todisk = TRUE)
options(warn = 1)
set.seed(123)

# ============================================================
# PATHS
# ============================================================

base_dir_plhiv <- DIR_PLHIV_FINE
base_dir_prev  <- DIR_PREV_FINE

out_dir        <- DIR_MC
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

clim_ref_path <- file.path(DIR_INDICES, "R95p", "ANNUAL_TOT", sprintf("R95pTOT_%d.tif", TARGET_YEAR))
gadm_path     <- FILE_GADM

# ============================================================
# SETTINGS
# ============================================================

year       <- as.character(TARGET_YEAR)
ages       <- c("15_19","20_24","25_29","30_34","35_39",
                "40_44","45_49","50_54","55_59")
K          <- length(ages)       # 9 age groups per sex
K_joint    <- 2 * K              # 18 total (9 female + 9 male)
z_crit     <- 1.96
eps        <- 1e-12
n_draws    <- as.integer(Sys.getenv("EWE_TEST_N_DRAWS", "1000"))
batch_size <- min(100L, n_draws)
datatype   <- "FLT4S"
heartbeat  <- TRUE

imputed_groups <- character(0)

# ============================================================
# Table S6 parameters and country-to-region mapping
# (Haeuser et al. 2022, BMC Medicine 20:488 - Supplementary Table S6)
# (Region assignments from Figure S12 of same paper)
# ============================================================

table_s6 <- data.frame(
  region     = c("Central", "Eastern", "Southern", "Western"),
  rho_age_Z2 = c(0.827258, 0.822233, 0.738138, 0.812046),
  rho_sex_Z2 = c(0.745520, 0.741752, 0.709504, 0.704386),
  var_Z2     = c(0.329295, 0.460780, 0.489873, 0.346227),
  rho_age_Z3 = c(0.838997, 0.991968, 0.797226, 0.884715),
  var_Z3     = c(0.094688, 0.912741, 0.051968, 0.182138),
  range_Z1   = c(0.083626, 0.028002, 0.046820, 0.051456),
  var_Z1     = c(0.227678, 0.310762, 0.123508, 0.224224),
  stringsAsFactors = FALSE
)

# Country -> region mapping (names match gadm_subset_countries.gpkg COUNTRY field)
country_to_region <- c(
  # Central SSA (6 countries)
  "Angola"                          = "Central",
  "Central African Republic"        = "Central",
  "Democratic Republic of the Congo"= "Central",
  "Equatorial Guinea"               = "Central",
  "Gabon"                           = "Central",
  "Republic of the Congo"           = "Central",

  # Eastern SSA (15 countries)
  "Burundi"       = "Eastern",
  "Djibouti"      = "Eastern",
  "Eritrea"       = "Eastern",
  "Ethiopia"      = "Eastern",
  "Kenya"         = "Eastern",
  "Madagascar"    = "Eastern",
  "Malawi"        = "Eastern",
  "Mozambique"    = "Eastern",
  "Rwanda"        = "Eastern",
  "Somalia"       = "Eastern",
  "South Sudan"   = "Eastern",
  "Sudan"         = "Eastern",
  "Tanzania"      = "Eastern",
  "Uganda"        = "Eastern",
  "Zambia"        = "Eastern",

  # Southern SSA (6 countries)
  "Botswana"      = "Southern",
  "Lesotho"       = "Southern",
  "Namibia"       = "Southern",
  "South Africa"  = "Southern",
  "Swaziland"     = "Southern",
  "Zimbabwe"      = "Southern",

  # Western SSA (16 countries)
  "Benin"                = "Western",
  "Burkina Faso"         = "Western",
  "Cameroon"             = "Western",
  "Chad"                 = "Western",
  "C\u00f4te d'Ivoire"  = "Western",
  "Gambia"               = "Western",
  "Ghana"                = "Western",
  "Guinea"               = "Western",
  "Guinea-Bissau"        = "Western",
  "Liberia"              = "Western",
  "Mali"                 = "Western",
  "Niger"                = "Western",
  "Nigeria"              = "Western",
  "Senegal"              = "Western",
  "Sierra Leone"         = "Western",
  "Togo"                 = "Western"
)

# ============================================================
# Correlation structures and variance weights
# ============================================================

message("Correlation structures and variance weights, Table S6")

# Per region:
#   Chol_Z2:     18x18 Cholesky of Z2 correlation (Kronecker)
#   Chol_Z3_age: 9x9 Cholesky of Z3 age correlation (AR(1))
#   w1, w2, w3:  variance fractions (sum to 1)

region_Chol_Z2    <- list()   # 18x18 Cholesky for Z2
region_Chol_Z3    <- list()   # 9x9 Cholesky for Z3 age component
region_w          <- list()   # variance fractions (w1, w2, w3)

for (r in seq_len(nrow(table_s6))) {
  reg <- table_s6$region[r]

  # Variance fractions
  v1   <- table_s6$var_Z1[r]
  v2   <- table_s6$var_Z2[r]
  v3   <- table_s6$var_Z3[r]
  vtot <- v1 + v2 + v3
  region_w[[reg]] <- c(w1 = v1/vtot, w2 = v2/vtot, w3 = v3/vtot)

  # Z2: Kronecker of sex(2x2) and age(KxK) AR(1)
  Sigma_age_Z2 <- toeplitz(table_s6$rho_age_Z2[r]^(0:(K-1)))
  Sigma_sex_Z2 <- matrix(c(1, table_s6$rho_sex_Z2[r],
                            table_s6$rho_sex_Z2[r], 1), 2, 2)
  Cor_Z2 <- kronecker(Sigma_sex_Z2, Sigma_age_Z2)
  region_Chol_Z2[[reg]] <- chol(Cor_Z2)

  # Z3: AR(1) across 9 age groups (applied independently to F and M)
  Cor_Z3_age <- toeplitz(table_s6$rho_age_Z3[r]^(0:(K-1)))
  region_Chol_Z3[[reg]] <- chol(Cor_Z3_age)

  w <- region_w[[reg]]
  message(sprintf("  %s: w1=%.1f%% (Z1 spatial), w2=%.1f%% (Z2 age-sex), w3=%.1f%% (Z3 country)",
                  reg, 100*w[1], 100*w[2], 100*w[3]))
  message(sprintf("         spatially shared: %.1f%%  |  Z2 cross-sex: %.3f  |  Z3 age-adj: %.3f",
                  100*(w[2]+w[3]), table_s6$rho_sex_Z2[r], table_s6$rho_age_Z3[r]))
}

# ============================================================
# Helper functions
# ============================================================

age_folder <- function(a) gsub("_", "-", a)

find_plhiv <- function(sex, bound, age) {
  subdir <- file.path(base_dir_plhiv, paste(sex, age_folder(age)),
                      toTitleCase(tolower(bound)))
  patt <- paste0("^", sex, "_", bound, "_", age, "_", year, "\\.tif$")
  f <- list.files(subdir, pattern = patt, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) {
    f <- list.files(base_dir_plhiv, pattern = patt, recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
  }
  if (length(f) == 0) NA_character_ else f[1]
}

find_prev <- function(sex, bound, age) {
  subdir <- file.path(base_dir_prev, paste(sex, age_folder(age)),
                      toTitleCase(tolower(bound)))
  patt <- paste0("HIV_PREVALENCE_", bound, "_", age, "_", sex, "_", year, ".*\\.TIF$")
  f <- list.files(subdir, pattern = patt, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) {
    f <- list.files(base_dir_prev, pattern = patt, recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
  }
  if (length(f) == 0 && toupper(bound) == "MEAN") {
    patt_imp <- paste0("HIV_PREVALENCE_MEAN_", age, "_", sex, "_", year, "_IMPUTED.*\\.TIF$")
    f <- list.files(base_dir_prev, pattern = patt_imp, recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
    if (length(f) == 0) {
      f <- list.files(out_dir, pattern = patt_imp, recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
    }
    if (length(f) > 0) {
      message(sprintf("  Imputed prevalence mean used for %s %s", sex, age))
    }
  }
  if (length(f) == 0) NA_character_ else f[1]
}

align_to_ref <- function(r, ref) {
  if (!same.crs(r, ref)) {
    r <- project(r, ref)
  } else if (!compareGeom(r, ref, stopOnError = FALSE)) {
    r <- resample(r, ref)
  }
  r
}

make_aligned_template <- function(climate_r, hiv_r) {
  if (is.character(climate_r)) climate_r <- rast(climate_r)
  if (is.character(hiv_r))     hiv_r     <- rast(hiv_r)
  tpl <- extend(climate_r, ext(hiv_r), snap = "out")
  tpl <- crop(tpl, ext(hiv_r), snap = "out")
  tpl[] <- 1:ncell(tpl)
  message(sprintf("  Template: %d cols x %d rows = %d cells @ %.4f deg",
                  ncol(tpl), nrow(tpl), ncell(tpl), res(tpl)[1]))
  tpl
}

verify_climate_alignment <- function(tpl, clim_ref, tolerance = 1e-9) {
  if (is.character(clim_ref)) clim_ref <- rast(clim_ref)
  res_ok  <- all(abs(res(tpl) - res(clim_ref)) < tolerance)
  x_off   <- (xmin(tpl) - xmin(clim_ref)) / res(clim_ref)[1]
  y_off   <- (ymin(tpl) - ymin(clim_ref)) / res(clim_ref)[2]
  grid_ok <- abs(x_off - round(x_off)) < tolerance &
    abs(y_off - round(y_off)) < tolerance
  crs_ok  <- same.crs(tpl, clim_ref)
  all_ok  <- res_ok && grid_ok && crs_ok
  if (all_ok) message("  Alignment checks passed")
  else        stop("  Alignment checks FAILED")
  invisible(all_ok)
}

build_overlap_W <- function(fine_ref, tpl, cache_rds = NULL) {
  if (!is.null(cache_rds) && file.exists(cache_rds)) {
    message("  Loading cached W: ", basename(cache_rds))
    obj <- readRDS(cache_rds)
    if (nrow(obj$W) == ncell(tpl)) return(obj)
    message("  Cached weights have the wrong dimensions, rebuilding")
  }
  if (is.character(fine_ref)) fine_ref <- rast(fine_ref)
  id <- rast(fine_ref); id[] <- 1:ncell(id)

  tpl_v <- deepcopy(tpl)
  tpl_v[] <- 1:ncell(tpl_v)

  polys <- st_as_sf(as.polygons(tpl_v, dissolve = FALSE))
  if (!same.crs(tpl_v, id)) {
    message("  CRS mismatch: reprojecting coarse polygons to fine grid CRS")
    polys <- st_transform(polys, crs(id))
  }

  message(sprintf("  Extracting overlaps: %d fine -> %d coarse...",
                  ncell(id), nrow(polys)))
  lst <- exact_extract(id, polys, include_xy = FALSE, progress = TRUE)

  nC  <- length(lst); nF <- ncell(id)
  lens <- vapply(lst, function(x) if (is.null(x)) 0L else nrow(x), integer(1))
  nnz  <- sum(lens)

  i <- integer(nnz); j <- integer(nnz); x <- numeric(nnz); pos <- 1L
  for (ci in seq_len(nC)) {
    df <- lst[[ci]]
    if (!is.null(df) && nrow(df) > 0) {
      n <- nrow(df); rng <- pos:(pos + n - 1L)
      i[rng] <- ci; j[rng] <- as.integer(df$value); x[rng] <- df$coverage_fraction
      pos <- pos + n
    }
  }
  if (pos <= nnz) { i <- i[1:(pos-1)]; j <- j[1:(pos-1)]; x <- x[1:(pos-1)] }

  W <- sparseMatrix(i = i, j = j, x = x, dims = c(nC, nF))
  obj <- list(W = W, row_sums = Matrix::rowSums(W), tpl = tpl)
  if (!is.null(cache_rds)) {
    saveRDS(obj, cache_rds)
    message("  Cached: ", basename(cache_rds))
  }
  obj
}

# ============================================================
# IMPUTATION
# ============================================================

impute_prevalence_mean <- function(sex, age) {
  message(sprintf("\n  *** IMPUTING prevalence MEAN for %s %s (lognormal + clamp) ***", sex, age))

  prev_L_path <- find_prev(sex, "LOWER", age)
  prev_U_path <- find_prev(sex, "UPPER", age)

  if (is.na(prev_L_path) || is.na(prev_U_path)) {
    stop(sprintf("Cannot impute %s %s: the lower or upper prevalence bound is also missing", sex, age))
  }

  r_L <- rast(prev_L_path)
  r_U <- rast(prev_U_path)

  if (!compareGeom(r_L, r_U, stopOnError = FALSE)) {
    r_U <- align_to_ref(r_U, r_L)
  }

  L_vals <- values(r_L, mat = FALSE)
  U_vals <- values(r_U, mat = FALSE)

  L_safe <- pmax(L_vals, 1e-8)
  U_safe <- pmax(U_vals, L_safe + 1e-8)

  mu_ln    <- (log(L_safe) + log(U_safe)) / 2
  sigma_ln <- (log(U_safe) - log(L_safe)) / (2 * 1.96)
  pred_mean <- pmin(exp(mu_ln + sigma_ln^2 / 2), 100)

  pred_mean[is.na(L_vals) | is.na(U_vals)] <- NA

  r_pred <- rast(r_L)
  values(r_pred) <- pred_mean

  imputed_path <- file.path(out_dir,
                            paste0("HIV_PREVALENCE_MEAN_", age, "_", sex, "_", year, "_IMPUTED.TIF"))
  writeRaster(r_pred, imputed_path, overwrite = TRUE,
              datatype = "FLT4S", gdal = c("COMPRESS=LZW", "PREDICTOR=3"))

  v <- pred_mean[!is.na(pred_mean)]
  n_clamped <- sum(v >= 99.99)
  message(sprintf("    Saved: %s", basename(imputed_path)))
  message(sprintf("    Range: %.6f - %.4f%%", min(v), max(v)))
  message(sprintf("    Median: %.4f%%  |  Non-NA pixels: %d", median(v), length(v)))
  if (n_clamped > 0) {
    message(sprintf("    Pixels clamped to 100%%: %d", n_clamped))
  }

  imputed_groups <<- c(imputed_groups, sprintf("%s_%s", sex, age))
  imputed_path
}

# ============================================================
# Check all files and impute if needed
# ============================================================

message("Checking input files and imputing missing means")

for (sex in c("FEMALES", "MALES")) {
  for (age in ages) {
    fL <- find_plhiv(sex, "LOWER", age)
    fM <- find_plhiv(sex, "MEAN",  age)
    fU <- find_plhiv(sex, "UPPER", age)
    if (anyNA(c(fL, fM, fU))) {
      stop(sprintf("MISSING PLHIV file for %s %s - cannot proceed", sex, age))
    }

    fPM <- find_prev(sex, "MEAN", age)
    if (is.na(fPM)) {
      fPL <- find_prev(sex, "LOWER", age)
      fPU <- find_prev(sex, "UPPER", age)
      if (!is.na(fPL) && !is.na(fPU)) {
        message(sprintf("  %s %s: Prevalence MEAN missing -> IMPUTING from L/U", sex, age))
        impute_prevalence_mean(sex, age)
      } else {
        stop(sprintf("MISSING prevalence MEAN, LOWER, and/or UPPER for %s %s", sex, age))
      }
    } else {
      message(sprintf("  %s %s: All files OK", sex, age))
    }
  }
}

if (length(imputed_groups) > 0) {
  message(sprintf("\n  IMPUTED GROUPS: %s", paste(imputed_groups, collapse = ", ")))
} else {
  message("\n  All prevalence means found - no imputation needed")
}
message("")

# ============================================================
# Step 1: template and overlap weights
# ============================================================

message("Step 1: template and overlap weights")

hiv_ref_path <- find_plhiv("FEMALES", "LOWER", ages[1])
if (is.na(hiv_ref_path)) stop("No PLHIV raster for the target year under DIR_PLHIV_FINE; run tests/00_preflight.R")
hiv_ref <- rast(hiv_ref_path)

tpl <- make_aligned_template(clim_ref_path, hiv_ref)
verify_climate_alignment(tpl, clim_ref_path)

nC <- ncell(tpl)
nF <- ncell(hiv_ref)

cache_W <- file.path(out_dir, sprintf("W_overlap_%dx%d_to_%dx%d.rds",
                                      ncol(hiv_ref), nrow(hiv_ref),
                                      ncol(tpl), nrow(tpl)))
Wobj <- build_overlap_W(hiv_ref, tpl, cache_rds = cache_W)
W    <- Wobj$W

fine_has_data_vec <- as.numeric(!is.na(values(rast(hiv_ref_path), mat = FALSE)))
coarse_has_data   <- as.numeric(W %*% fine_has_data_vec) > 0  # preliminary; recomputed from W_data later
n_data_cells      <- sum(coarse_has_data)

message(sprintf("  Fine cells: %d | Coarse cells: %d", nF, nC))
message(sprintf("  Coarse cells with data: %d / %d (%.1f%%)",
                n_data_cells, nC, 100 * n_data_cells / nC))

# ============================================================
# Step 2: fine-resolution shifted lognormal parameters and population
#         X = gamma + exp(mu + sigmaZ)   three parameters fitted to L, M, U exactly.
#         sigma from ratio R=(M-L)/(U-L) via lookup, then mu and gamma analytically.
# ============================================================

message("Step 2: shifted lognormal parameters and aggregated population")

# Pre-compute lookup table for R(sigma) -> sigma inversion
# R(sigma) = (exp(sigma2/2) - exp(-1.96sigma)) / (exp(1.96sigma) - exp(-1.96sigma))
# R is monotonically decreasing in sigma.
sg_grid <- seq(0.001, 2.0, length.out = 20000)
R_lookup <- (exp(sg_grid^2 / 2) - exp(-1.96 * sg_grid)) /
            (exp(1.96 * sg_grid) - exp(-1.96 * sg_grid))
sg_grid_rev <- rev(sg_grid)
R_lookup_rev <- rev(R_lookup)

message(sprintf("  R(sigma) lookup table: %d points, sigma in [%.3f, %.3f], R in [%.4f, %.4f]",
                length(sg_grid), min(sg_grid), max(sg_grid), min(R_lookup), max(R_lookup)))

compute_fine_params_both_sexes <- function() {

  fine_mu_all    <- vector("list", K_joint)
  fine_sg_all    <- vector("list", K_joint)
  fine_gm_all    <- vector("list", K_joint)   # gamma (shift parameter)
  fine_zmask_all <- vector("list", K_joint)
  coarse_POP_all <- vector("list", K_joint)
  PLHIV_mean_all <- vector("list", K_joint)

  joint_names <- c(paste0("F_", ages), paste0("M_", ages))

  for (sx_idx in 1:2) {
    sex <- c("FEMALES", "MALES")[sx_idx]
    offset <- (sx_idx - 1) * K

    message(sprintf("\n--- %s (columns %d-%d of 18) ---", sex, offset + 1, offset + K))

    message("  Loading fine-resolution PLHIV L, M, U...")
    Lp <- lapply(ages, function(a) {
      f <- find_plhiv(sex, "LOWER", a)
      if (is.na(f)) stop(sprintf("Missing PLHIV LOWER for %s %s", sex, a))
      values(rast(f), mat = FALSE)
    })
    Up <- lapply(ages, function(a) {
      f <- find_plhiv(sex, "UPPER", a)
      if (is.na(f)) stop(sprintf("Missing PLHIV UPPER for %s %s", sex, a))
      values(rast(f), mat = FALSE)
    })
    Mp <- lapply(ages, function(a) {
      f <- find_plhiv(sex, "MEAN", a)
      if (is.na(f)) stop(sprintf("Missing PLHIV MEAN for %s %s", sex, a))
      values(rast(f), mat = FALSE)
    })

    message("  Loading fine-resolution prevalence means...")
    Mv <- lapply(ages, function(a) {
      f <- find_prev(sex, "MEAN", a)
      if (is.na(f)) stop(sprintf("Missing PREVALENCE MEAN for %s %s", sex, a))
      values(rast(f), mat = FALSE)
    })

    message("  Fitting shifted lognormal parameters at fine resolution...")
    for (i in seq_len(K)) {
      idx <- offset + i

      L <- pmax(Lp[[i]], eps)
      U <- pmax(Up[[i]], eps)
      U <- pmax(U, L + eps)
      M <- pmax(Mp[[i]], eps)
      M <- pmax(M, L + eps)
      M <- pmin(M, U - eps)

      R_data <- (M - L) / (U - L)
      sg <- approx(R_lookup_rev, sg_grid_rev, xout = R_data, rule = 2)$y

      exp_mu <- (U - L) / (exp(1.96 * sg) - exp(-1.96 * sg))
      exp_mu <- pmax(exp_mu, eps)
      mu <- log(exp_mu)
      gm <- L - exp(mu - 1.96 * sg)

      # Cells where L, M and U are effectively equal carry no uncertainty to
      # propagate, so the draw is fixed at the mean.
      degenerate <- !is.finite(sg) | !is.finite(mu) | !is.finite(gm) | sg < 1e-8
      n_degen <- sum(degenerate & !is.na(Lp[[i]]))
      if (any(degenerate, na.rm = TRUE)) {
        # For degenerate cells: gamma = 0, standard lognormal (no shift needed)
        # These cells have near-zero uncertainty (L ~= M ~= U), so the
        # distribution choice barely matters   draws cluster tightly around M.
        sg[degenerate] <- pmax((log(U[degenerate]) - log(L[degenerate])) / (2 * z_crit), 1e-6)
        mu[degenerate] <- log(M[degenerate]) - sg[degenerate]^2 / 2
        gm[degenerate] <- 0
      }

      # Treat near-deterministic cells where published L ~= M ~= U
      near_determ <- (!is.na(Lp[[i]])) &
                     ((Up[[i]] - Lp[[i]]) < 1e-6 |
                      (abs(Up[[i]] - Mp[[i]]) < 1e-6 & abs(Mp[[i]] - Lp[[i]]) < 1e-6))
      n_determ <- sum(near_determ & !degenerate)
      if (any(near_determ)) {
        sg[near_determ] <- 1e-6
        mu[near_determ] <- log(M[near_determ])
        gm[near_determ] <- 0
      }

      fine_mu_all[[idx]] <- mu
      fine_sg_all[[idx]] <- sg
      fine_gm_all[[idx]] <- gm

      fine_zmask_all[[idx]] <- Up[[i]] <= 0

      fine_mu_all[[idx]][is.na(fine_mu_all[[idx]])] <- 0
      fine_sg_all[[idx]][is.na(fine_sg_all[[idx]])] <- 0
      fine_gm_all[[idx]][is.na(fine_gm_all[[idx]])] <- 0

      n_zero <- sum(fine_zmask_all[[idx]] & !is.na(Lp[[i]]))
      n_eps  <- sum(Lp[[i]] == 0 & !fine_zmask_all[[idx]], na.rm = TRUE)
      if (n_zero > 0 || n_eps > 0 || n_degen > 0 || n_determ > 0) {
        message(sprintf("    %s: %d zero-masked, %d L=0, %d degenerate, %d deterministic",
                        joint_names[idx], n_zero, n_eps, n_degen, n_determ))
      }

      # -- Reconstruction diagnostic: verify (gamma,mu,sigma) recover published L/M/U --
      # Skip NA and degenerate/deterministic cells
      ok <- !is.na(Lp[[i]]) & !degenerate & !near_determ & !fine_zmask_all[[idx]]
      if (sum(ok) > 0) {
        Lhat <- gm[ok] + exp(mu[ok] - 1.96 * sg[ok])
        Mhat <- gm[ok] + exp(mu[ok] + sg[ok]^2 / 2)
        Uhat <- gm[ok] + exp(mu[ok] + 1.96 * sg[ok])
        err_L <- max(abs(Lhat / L[ok] - 1))
        err_M <- max(abs(Mhat / M[ok] - 1))
        err_U <- max(abs(Uhat / U[ok] - 1))
        if (max(err_L, err_M, err_U) > 0.001) {
          message(sprintf("    %s: reconstruction max error L=%.4f%% M=%.4f%% U=%.4f%%",
                          joint_names[idx], 100*err_L, 100*err_M, 100*err_U))
        }
      }

      PLHIV_mean_all[[idx]] <- Mp[[i]]
    }

    message("  Computing and aggregating population per age group...")
    for (i in seq_len(K)) {
      idx <- offset + i

      prev_prop <- Mv[[i]] / 100
      prev_prop[is.na(prev_prop) | prev_prop <= 0] <- NA

      fine_pop <- Mp[[i]] / prev_prop
      fine_pop[!is.finite(fine_pop)] <- 0
      fine_pop[is.na(fine_pop)]      <- 0

      coarse_POP_all[[idx]] <- as.numeric(W %*% fine_pop)

      grp_label <- sprintf("%s_%s", sex, ages[i])
      if (grp_label %in% imputed_groups) {
        message(sprintf("    %s: POP uses IMPUTED prevalence mean", joint_names[idx]))
      }
    }
  }

  names(fine_mu_all)    <- joint_names
  names(fine_sg_all)    <- joint_names
  names(fine_gm_all)    <- joint_names
  names(fine_zmask_all) <- joint_names
  names(coarse_POP_all) <- joint_names
  names(PLHIV_mean_all) <- joint_names

  # Data mask: a cell must carry a value in all 18 age-sex groups.

  fine_data_mask <- rep(TRUE, nF)
  for (idx in seq_len(K_joint)) {
    fine_data_mask <- fine_data_mask & !is.na(PLHIV_mean_all[[idx]])
  }

  data_idx <- which(fine_data_mask)
  nD <- length(data_idx)
  message(sprintf("\n  Fine cells with data (all 18 groups): %d / %d (%.1f%%)",
                  nD, nF, 100 * nD / nF))

  fine_mu_d    <- lapply(fine_mu_all,    function(v) v[data_idx])
  fine_sg_d    <- lapply(fine_sg_all,    function(v) v[data_idx])
  fine_gm_d    <- lapply(fine_gm_all,    function(v) v[data_idx])
  fine_zmask_d <- lapply(fine_zmask_all, function(v) v[data_idx])
  names(fine_mu_d)    <- joint_names
  names(fine_sg_d)    <- joint_names
  names(fine_gm_d)    <- joint_names
  names(fine_zmask_d) <- joint_names

  W_data <- W[, data_idx, drop = FALSE]
  message(sprintf("  W subsetted: %d x %d -> %d x %d",
                  nrow(W), ncol(W), nrow(W_data), ncol(W_data)))

  coarse_POP_F_total <- Reduce("+", coarse_POP_all[1:K])
  coarse_POP_M_total <- Reduce("+", coarse_POP_all[(K+1):K_joint])
  coarse_POP_ALL_total <- coarse_POP_F_total + coarse_POP_M_total

  message(sprintf("  POP FEMALES total: %.0f", sum(coarse_POP_F_total[coarse_has_data])))
  message(sprintf("  POP MALES total:   %.0f", sum(coarse_POP_M_total[coarse_has_data])))
  message(sprintf("  POP ALL total:     %.0f", sum(coarse_POP_ALL_total[coarse_has_data])))

  list(
    fine_mu_d    = fine_mu_d,
    fine_sg_d    = fine_sg_d,
    fine_gm_d    = fine_gm_d,
    fine_zmask_d = fine_zmask_d,
    nD           = nD,
    data_idx     = data_idx,
    W_data       = W_data,
    POP_F_total  = coarse_POP_F_total,
    POP_M_total  = coarse_POP_M_total,
    POP_ALL_total = coarse_POP_ALL_total
  )
}

params <- compute_fine_params_both_sexes()

# Recompute coarse_has_data from W_data (not from single reference raster)
coarse_has_data <- rowSums(params$W_data) > 0
n_data_cells    <- sum(coarse_has_data)
message(sprintf("  coarse_has_data (from W_data): %d / %d cells", n_data_cells, nC))

# ============================================================
# Step 3: assign fine data cells to countries and regions
#
# The spatial join also builds:
#   - country_idx: list mapping country name -> fine cell indices
#   - cell_region_name: region of each fine cell (for w1/w2/w3 lookup)
#   - country_region: which region each country belongs to
# These are needed for the three-component Z generation in Step 4.
# ============================================================

message("Step 3: assign fine cells to countries and regions VIA GADM")

# Read GADM boundaries
countries_sf <- st_read(gadm_path, quiet = TRUE)
countries_sf <- st_make_valid(countries_sf)

# Get XY coordinates of data cells in WGS84 lon/lat (degrees)
# The Z1 range in Table S6 is in degrees, so coordinates are kept in degrees.
# Always transform explicitly to EPSG:4326   no guessing from CRS string.
xy_native <- xyFromCell(hiv_ref, params$data_idx)
hiv_crs <- crs(hiv_ref)
tmp_pts <- st_as_sf(data.frame(x = xy_native[, 1], y = xy_native[, 2]),
                    coords = c("x", "y"), crs = hiv_crs)
tmp_pts_4326 <- st_transform(tmp_pts, "EPSG:4326")
xy <- st_coordinates(tmp_pts_4326)  # guaranteed degrees (lon, lat)

pts <- st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                coords = c("x", "y"), crs = "EPSG:4326")

# Reproject if needed
if (!identical(st_crs(pts), st_crs(countries_sf))) {
  message("  Reprojecting points to GADM CRS...")
  pts <- st_transform(pts, st_crs(countries_sf))
}

# Spatial join - disable S2 to avoid invalid polygon errors
message(sprintf("  Joining %d data cells to country polygons...", params$nD))
s2_was_on <- sf::sf_use_s2()
sf::sf_use_s2(FALSE)
joined <- st_join(pts, countries_sf[, "COUNTRY"], join = st_intersects)
sf::sf_use_s2(s2_was_on)
cell_countries <- joined$COUNTRY

# Map country names to regions
cell_regions <- country_to_region[cell_countries]

# Handle unmatched cells - two distinct cases:
#   A. a country name that is absent from country_to_region, which stops the step
#   B. no country at all, for ocean and borders, resolved to the nearest country

# Case A: named countries not in mapping
named_but_unmatched <- !is.na(cell_countries) & is.na(cell_regions)
if (any(named_but_unmatched)) {
  bad_countries <- unique(cell_countries[named_but_unmatched])
  message("\n  Country names absent from the country_to_region mapping:")
  for (bc in bad_countries) {
    n_bc <- sum(cell_countries == bc, na.rm = TRUE)
    message(sprintf("    '%s' (%d cells)", bc, n_bc))
  }
  stop("Add these countries to the country_to_region mapping and run again.")
}

# Case B: no country from GADM (typically <50 border/ocean)
no_country <- is.na(cell_countries)
n_no_country <- sum(no_country)
if (n_no_country > 0) {
  message(sprintf("  %d cells had no GADM country (border/ocean slivers)", n_no_country))

  # Assign each to nearest country via spatial nearest-neighbour
  orphan_pts <- pts[no_country, ]
  matched_pts <- pts[!no_country, ]
  nearest_idx <- st_nearest_feature(orphan_pts, matched_pts)
  cell_countries[no_country] <- cell_countries[which(!no_country)][nearest_idx]
  cell_regions[no_country] <- country_to_region[cell_countries[no_country]]

  message(sprintf("  Assigned to nearest country: %s",
                  paste(unique(cell_countries[no_country]), collapse = ", ")))

  if (any(is.na(cell_regions))) {
    stop("Still have unmatched cells after nearest-neighbour. Check GADM/mapping.")
  }
}

# Build index vectors per region
region_idx <- list()
for (reg in table_s6$region) {
  region_idx[[reg]] <- which(cell_regions == reg)
}

# Build index vectors per country
# All cells have valid country names at this point 
unique_countries <- sort(unique(cell_countries))
country_idx <- list()
for (cname in unique_countries) {
  country_idx[[cname]] <- which(cell_countries == cname)
}

# Map each country to its region (for Cholesky lookup)
country_region <- country_to_region[unique_countries]

# Report
region_table <- table(cell_regions)
message("\n  Region assignment:")
for (reg in names(region_table)) {
  n_countries_in_reg <- sum(country_region == reg)
  message(sprintf("    %s: %d cells (%.1f%%), %d countries",
                  reg, region_table[reg],
                  100 * region_table[reg] / params$nD,
                  n_countries_in_reg))
}
message(sprintf("\n  Total countries with data cells: %d", length(unique_countries)))

# -- Z1 SPATIAL BLOCKS --
# Coordinates are discretised onto a grid whose spacing is the Matern range
# from Table S6, and fine cells in the same block share one Z1 value. Cells
# within a range are then highly correlated and cells beyond it independent.
#

message("\n  Building Z1 spatial blocks (Matern range blocks)...")

# Get region-specific ranges
region_range <- setNames(table_s6$range_Z1, table_s6$region)

# For each fine cell, compute which Z1 block it belongs to
# Block ID = unique combination of (region, floor(x/range), floor(y/range))
z1_block_id <- integer(params$nD)   # maps fine cell index -> block ID
block_counter <- 0L

z1_n_blocks <- list()  # per region: how many blocks

for (reg in table_s6$region) {
  idx <- region_idx[[reg]]
  if (length(idx) == 0) next

  rng <- region_range[reg]
  xy_reg <- xy[idx, , drop = FALSE]  # xy from Step 3 (all data cells)

  # Discretise to blocks
  bx <- floor(xy_reg[, 1] / rng)
  by <- floor(xy_reg[, 2] / rng)

  # Map (bx, by) pairs to unique block IDs within this region
  block_key <- paste(bx, by, sep = "_")
  unique_keys <- unique(block_key)
  n_blocks_reg <- length(unique_keys)
  key_to_local <- match(block_key, unique_keys)

  # Offset by global counter so block IDs are unique across regions
  z1_block_id[idx] <- key_to_local + block_counter
  block_counter <- block_counter + n_blocks_reg
  z1_n_blocks[[reg]] <- n_blocks_reg

  message(sprintf("    %s: range=%.4f deg, %d blocks for %d cells (avg %.1f cells/block)",
                  reg, rng, n_blocks_reg, length(idx), length(idx) / max(n_blocks_reg, 1)))
}

n_z1_blocks <- block_counter
message(sprintf("  Total Z1 blocks: %d (vs %d fine cells = %.1fx compression)",
                n_z1_blocks, params$nD, params$nD / max(n_z1_blocks, 1)))

# ============================================================
# Step 4: Monte Carlo sampling, three-component decomposition
#
# Per draw:
#   1. Z1  one value per spatial block, shared across all 18 age-sex groups
#   2. Z2  one 18-dimensional draw per region, constant within the region
#   3. Z3  one 9-dimensional draw per country per sex
#
#   Z_total[p, k] = sqrt(w1) Z1[p] + sqrt(w2) Z2_r[k] + sqrt(w3) Z3_c[k]
#   PLHIV[p, k]   = gamma[p,k] + exp(mu[p,k] + sigma[p,k] Z_total[p,k])
#
# Z2 and Z3 are shared spatially and so survive aggregation of fine cells.
#
# ============================================================

message("Step 4: Monte Carlo sampling")

fine_mu_d    <- params$fine_mu_d
fine_sg_d    <- params$fine_sg_d
fine_gm_d    <- params$fine_gm_d
fine_zmask_d <- params$fine_zmask_d
nD           <- params$nD
W_data       <- params$W_data
POP_F_total  <- params$POP_F_total
POP_M_total  <- params$POP_M_total
POP_ALL_total <- params$POP_ALL_total

n_batches <- ceiling(n_draws / batch_size)

# Pre-compute sqrt(w) per fine cell for combining
# Each fine cell inherits its region's variance fractions
sqrt_w1 <- numeric(nD)
sqrt_w2 <- numeric(nD)
sqrt_w3 <- numeric(nD)
for (reg in names(region_idx)) {
  idx <- region_idx[[reg]]
  w <- region_w[[reg]]
  sqrt_w1[idx] <- sqrt(w["w1"])
  sqrt_w2[idx] <- sqrt(w["w2"])
  sqrt_w3[idx] <- sqrt(w["w3"])
}

message(sprintf("  Fine data cells: %d | Coarse cells: %d", nD, nC))
message(sprintf("  Countries: %d | Regions: %d", length(unique_countries), length(region_idx)))
message(sprintf("  Draws: %d | Batches: %d", n_draws, n_batches))
message("  Strategy: Z1 (per spatial block) + Z2 (per region) + Z3 (per country)")
message(sprintf("  Z1 blocks: %d (vs %d fine cells)", n_z1_blocks, nD))

for (b in seq_len(n_batches)) {
  nb <- min(batch_size, n_draws - (b - 1) * batch_size)

  batch_file_F   <- file.path(out_dir, sprintf("FEMALES_15_59_%s_b%03d.tif", year, b))
  batch_file_M   <- file.path(out_dir, sprintf("MALES_15_59_%s_b%03d.tif", year, b))
  batch_file_ALL <- file.path(out_dir, sprintf("ALL_15_59_%s_b%03d.tif", year, b))

  if (file.exists(batch_file_F) && file.exists(batch_file_M) && file.exists(batch_file_ALL)) {
    message(sprintf("  [Batch %03d/%03d] All exist -> skipping", b, n_batches))
    next
  }

  message(sprintf("  [Batch %03d/%03d] Processing %d draws...", b, n_batches, nb))
  t0 <- Sys.time()

  out_layers_F   <- vector("list", nb)
  out_layers_M   <- vector("list", nb)
  out_layers_ALL <- vector("list", nb)

  for (k in seq_len(nb)) {
    draw_id <- (b - 1) * batch_size + k

    # -- 1. Z1: one scalar per SPATIAL BLOCK (shared across all 18 groups) --
    # Fine cells in the same block (within one Matern range) share the same Z1.
    # This preserves spatial correlation at the correct scale per region.
    z1_blocks <- rnorm(n_z1_blocks)          # one draw per block
    Z1 <- z1_blocks[z1_block_id]             # map to fine cells

    # -- 2. Z2: one 18-dim draw per region --
    # Z2_expanded[p, k] = Z2_r[k] for each fine cell p in region r
    Z2_expanded <- matrix(0, nrow = nD, ncol = K_joint)
    for (reg in names(region_idx)) {
      idx <- region_idx[[reg]]
      if (length(idx) == 0) next
      z2_r <- as.numeric(rnorm(K_joint) %*% region_Chol_Z2[[reg]])  # 18-dim correlated
      Z2_expanded[idx, ] <- matrix(z2_r, nrow = length(idx), ncol = K_joint, byrow = TRUE)
    }

    # -- 3. Z3: one 9-dim draw per country, SEPARATE for F and M --
    # Z3_expanded[p, k] = Z3_c[k] for each fine cell p in country c
    # Z3 is assumed independent between sexes
    Z3_expanded <- matrix(0, nrow = nD, ncol = K_joint)
    for (cname in unique_countries) {
      cidx <- country_idx[[cname]]
      if (length(cidx) == 0) next
      reg <- country_region[[cname]]

      # Independent draws for females (cols 1-9) and males (cols 10-18)
      z3_F <- as.numeric(rnorm(K) %*% region_Chol_Z3[[reg]])   # 9-dim correlated
      z3_M <- as.numeric(rnorm(K) %*% region_Chol_Z3[[reg]])   # independent 9-dim
      z3_18 <- c(z3_F, z3_M)
      Z3_expanded[cidx, ] <- matrix(z3_18, nrow = length(cidx), ncol = K_joint, byrow = TRUE)
    }

    # -- 4. Combine: Z_total = sqrt(w1)*Z1 + sqrt(w2)*Z2 + sqrt(w3)*Z3 --
    # Z1 is scalar per cell -> broadcasts to all 18 columns
    Z_total <- sweep(Z2_expanded, 1, sqrt_w2, "*") +
               sweep(Z3_expanded, 1, sqrt_w3, "*")
    Z_total <- Z_total + (sqrt_w1 * Z1)  # Z1 adds same value to all 18 columns

    # -- 5. Sample FEMALE age groups (columns 1-9), sum --
    #    SHIFTED LOGNORMAL: PLHIV = gamma + exp(mu + sigma  x  Z_total)
    fine_plhiv_F <- numeric(nD)
    for (i in 1:K) {
      plhiv_i <- fine_gm_d[[i]] + exp(fine_mu_d[[i]] + fine_sg_d[[i]] * Z_total[, i])
      plhiv_i[fine_zmask_d[[i]]] <- 0
      plhiv_i[plhiv_i < 0] <- 0
      plhiv_i[!is.finite(plhiv_i)] <- 0
      fine_plhiv_F <- fine_plhiv_F + plhiv_i
    }

    # -- 6. Sample MALE age groups (columns 10-18), sum --
    #    SHIFTED LOGNORMAL: PLHIV = gamma + exp(mu + sigma  x  Z_total)
    fine_plhiv_M <- numeric(nD)
    for (i in 1:K) {
      j <- K + i
      plhiv_i <- fine_gm_d[[j]] + exp(fine_mu_d[[j]] + fine_sg_d[[j]] * Z_total[, j])
      plhiv_i[fine_zmask_d[[j]]] <- 0
      plhiv_i[plhiv_i < 0] <- 0
      plhiv_i[!is.finite(plhiv_i)] <- 0
      fine_plhiv_M <- fine_plhiv_M + plhiv_i
    }

    # -- 7. Combined = female + male --
    fine_plhiv_ALL <- fine_plhiv_F + fine_plhiv_M

    # -- 8. Aggregate to coarse via W_data --
    coarse_plhiv_F   <- as.numeric(W_data %*% fine_plhiv_F)
    coarse_plhiv_M   <- as.numeric(W_data %*% fine_plhiv_M)
    coarse_plhiv_ALL <- as.numeric(W_data %*% fine_plhiv_ALL)

    coarse_plhiv_F[!coarse_has_data]   <- NA
    coarse_plhiv_M[!coarse_has_data]   <- NA
    coarse_plhiv_ALL[!coarse_has_data] <- NA

    pop_F   <- POP_F_total;   pop_F[!coarse_has_data]   <- NA
    pop_M   <- POP_M_total;   pop_M[!coarse_has_data]   <- NA
    pop_ALL <- POP_ALL_total; pop_ALL[!coarse_has_data] <- NA

    # -- 9. Derive prevalence --
    derive_prev <- function(plhiv, pop) {
      prev <- ifelse(
        is.na(plhiv) | is.na(pop), NA,
        ifelse(pop == 0 | plhiv == 0, 0, 100 * plhiv / pop)
      )
      pmin(pmax(prev, 0, na.rm = FALSE), 100)
    }

    prev_F   <- derive_prev(coarse_plhiv_F,   pop_F)
    prev_M   <- derive_prev(coarse_plhiv_M,   pop_M)
    prev_ALL <- derive_prev(coarse_plhiv_ALL, pop_ALL)

    # -- Build rasters --
    build_draw_stack <- function(plhiv, pop, prev, draw_id) {
      r_pl  <- rast(tpl); values(r_pl)  <- plhiv
      r_pop <- rast(tpl); values(r_pop) <- pop
      r_pr  <- rast(tpl); values(r_pr)  <- prev
      names(r_pl)  <- sprintf("PLHIV_d%04d", draw_id)
      names(r_pop) <- sprintf("POP_d%04d", draw_id)
      names(r_pr)  <- sprintf("PREVpct_d%04d", draw_id)
      c(r_pl, r_pop, r_pr)
    }

    out_layers_F[[k]]   <- build_draw_stack(coarse_plhiv_F,   pop_F,   prev_F,   draw_id)
    out_layers_M[[k]]   <- build_draw_stack(coarse_plhiv_M,   pop_M,   prev_M,   draw_id)
    out_layers_ALL[[k]] <- build_draw_stack(coarse_plhiv_ALL, pop_ALL, prev_ALL, draw_id)

    if (heartbeat && (k %% 10 == 0 || k == nb)) {
      message(sprintf("    draw %d/%d done", k, nb))
    }
  }

  write_batch <- function(layers, path) {
    stk <- do.call(c, layers)
    writeRaster(stk, path, overwrite = TRUE, datatype = datatype,
                gdal = c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6"))
  }

  write_batch(out_layers_F,   batch_file_F)
  write_batch(out_layers_M,   batch_file_M)
  write_batch(out_layers_ALL, batch_file_ALL)

  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("  Batch %03d written (%.1f min): F + M + ALL", b, dt))

  rm(out_layers_F, out_layers_M, out_layers_ALL, Z_total, Z2_expanded, Z3_expanded); gc(verbose = FALSE)
}

message("\n  All batches complete\n")

# ============================================================
# Step 5: stitch batches and compute uncertainty intervals
#
# ============================================================

stitch_batches <- function(group_label, n_draws_target = 1000) {
  message(sprintf("\n=== STITCHING: %s ===", group_label))

  patt <- sprintf("^%s_15_59_%s_b[0-9]{3}\\.tif$", group_label, year)
  files <- sort(list.files(out_dir, pattern = patt, full.names = TRUE))
  if (!length(files)) stop("Monte Carlo batches are missing, the sampling stage did not complete for ", group_label)

  message(sprintf("  Loading %d batches...", length(files)))
  all_draws <- do.call(c, lapply(files, rast))

  idx_pl   <- grep("^PLHIV_d\\d+$", names(all_draws))
  ids_pl   <- as.integer(sub("^PLHIV_d(\\d+)$", "\\1", names(all_draws)[idx_pl]))
  keep     <- !duplicated(ids_pl)
  kept_ids <- ids_pl[keep]
  if (length(kept_ids) > n_draws_target) {
    kept_ids <- kept_ids[seq_len(n_draws_target)]
  }

  all_names <- names(all_draws)
  idx_keep <- integer(0)
  for (id in kept_ids) {
    idx_keep <- c(idx_keep,
                  match(sprintf("PLHIV_d%04d", id),   all_names),
                  match(sprintf("POP_d%04d", id),      all_names),
                  match(sprintf("PREVpct_d%04d", id),  all_names))
  }
  idx_keep <- idx_keep[!is.na(idx_keep)]

  dedup <- all_draws[[idx_keep]]

  full_path <- file.path(out_dir,
                         sprintf("%s_15_59_%s_draws_0p25.tif", group_label, year))
  writeRaster(dedup, full_path, overwrite = TRUE, datatype = datatype,
              gdal = c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6"))
  message("  Written: ", basename(full_path))
  message(sprintf("  Total draws: %d | Bands: %d", length(kept_ids), nlyr(dedup)))

  # Paper uses mean of 1,000 draws as point estimate 
  # with 2.5th and 97.5th percentiles for 95% UI
  qfun <- function(v) {
    c(lower = quantile(v, 0.025, na.rm = TRUE),
      mean  = mean(v, na.rm = TRUE),
      upper = quantile(v, 0.975, na.rm = TRUE))
  }

  idx_plhiv <- grep("^PLHIV_d",   names(dedup))
  idx_pop   <- grep("^POP_d",     names(dedup))
  idx_prev  <- grep("^PREVpct_d", names(dedup))

  ui_plhiv <- app(dedup[[idx_plhiv]], qfun)
  names(ui_plhiv) <- c("PLHIV_lower", "PLHIV_mean", "PLHIV_upper")

  ui_pop <- app(dedup[[idx_pop]], qfun)
  names(ui_pop) <- c("POP_lower", "POP_mean", "POP_upper")

  ui_prev <- app(dedup[[idx_prev]], qfun)
  names(ui_prev) <- c("PREV_lower", "PREV_mean", "PREV_upper")

  ui_stack <- c(ui_plhiv, ui_pop, ui_prev)
  ui_path <- file.path(out_dir,
                       sprintf("%s_15_59_%s_UI_0p25.tif", group_label, year))
  writeRaster(ui_stack, ui_path, overwrite = TRUE, datatype = datatype,
              gdal = c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6"))
  message("  UI written: ", basename(ui_path))

  pl_mean <- values(ui_plhiv[[2]], mat = FALSE)
  pl_mean <- pl_mean[!is.na(pl_mean)]
  message(sprintf("  PLHIV mean: total = %.0f, cell range [%.1f, %.1f]",
                  sum(pl_mean), min(pl_mean), max(pl_mean)))
}

stitch_batches("FEMALES")
stitch_batches("MALES")
stitch_batches("ALL")

# ============================================================
# Step 6: country layer
# ============================================================

message("Step 6: country layer")

ci     <- compute_country_raster(tpl, gadm_path)
cnames <- sort(unique(na.omit(ci$assignment)))
cmap   <- setNames(seq_along(cnames), cnames)
r_country <- tpl; values(r_country) <- cmap[ci$assignment]
names(r_country) <- "country"
levels(r_country) <- data.frame(value = seq_along(cnames), country = cnames)

for (f_name in c(sprintf("FEMALES_15_59_%s_draws_0p25.tif", year),
                 sprintf("MALES_15_59_%s_draws_0p25.tif", year),
                 sprintf("ALL_15_59_%s_draws_0p25.tif", year))) {
  f_path <- file.path(out_dir, f_name)
  if (file.exists(f_path)) {
    r <- rast(f_path)
    r_out <- c(r, r_country)
    out_f <- sub("\\.tif$", "_withCountry.tif", f_path)
    writeRaster(r_out, out_f, overwrite = TRUE, datatype = datatype,
                gdal = c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6"))
    message("  Written: ", basename(out_f))
  }
}

# ============================================================
# Aggregate uncertainty intervals
# ============================================================

for (grp in c("FEMALES", "MALES", "ALL")) {
  d <- rast(file.path(out_dir, sprintf("%s_15_59_%s_draws_0p25.tif", grp, year)))
  tot <- global(d[[grep("^PLHIV_d", names(d))]], "sum", na.rm = TRUE)$sum
  ci  <- quantile(tot, c(0.025, 0.975))
  message(sprintf("  %-8s %12.0f  [%.0f, %.0f]  width %.1f%%",
                  grp, mean(tot), ci[1], ci[2], 100 * (ci[2] - ci[1]) / mean(tot)))
}

message("Done. Next: step 13.")
