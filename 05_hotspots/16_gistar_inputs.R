# ============================================================
# Gi* input stack
# ------------------------------------------------------------
# Assembles the prevalence draws, the country layer and the nine indicators
# into one stack, applies the sign convention so that higher values denote
# worse conditions, and applies a cube root.
#
# Requires: steps 02 to 12
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))


# Stage 1 writes an augmented stack: the prevalence draws, the country layer
# and the nine indicator surfaces on one grid.
#
#   rainfall           R95pTOT, R99pTOT
#   heatwaves          TX90p_HWMF, TX95p_HWMF
#   drought at <= -1.0 SPI3, SPEI3, SRI3, SMA3
#   drought at <= -1.5 SPI3
#
# Stage 2 writes the stack that step 18 reads. Drought severity sums are
# negated so that higher values denote worse conditions across every
# indicator, and a sign-preserving cube root, sign(x) * |x|^(1/3), is then
# applied to every layer except the country layer. No standardisation is
# applied, since Gi* standardises internally.
#
# The transformation is the one selected in step 17.

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

# ============================================================
# CONFIG
# ============================================================

config <- list(
  year = TARGET_YEAR,

  # Input stacks from step 12, with the country layer appended
  stacks = list(
    ALL = file.path(DIR_MC, sprintf("ALL_15_59_%d_draws_0p25_withCountry.tif", TARGET_YEAR))
  ),


  # -- Rainfall --
  rainfall = list(
    R95pTOT = file.path(DIR_INDICES, sprintf("R95p/ANNUAL_TOT/R95pTOT_%d.tif", TARGET_YEAR)),
    R99pTOT = file.path(DIR_INDICES, sprintf("R99p/ANNUAL_TOT/R99pTOT_%d.tif", TARGET_YEAR))
  ),

  # -- Heatwaves --
  heatwaves = list(
    TX90p_HWMF = file.path(DIR_INDICES, sprintf("Heatwaves/TX90p/ANN/HWMF/TX90p_ANN_HWMF_%d.tif", TARGET_YEAR)),
    TX95p_HWMF = file.path(DIR_INDICES, sprintf("Heatwaves/TX95p/ANN/HWMF/TX95p_ANN_HWMF_%d.tif", TARGET_YEAR))
  ),

  # -- Drought severity sums --
  drought_leq1p0 = list(
    SPI3  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
    SPEI3 = file.path(DIR_INDICES, sprintf("SPEI/SPEI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
    SRI3  = file.path(DIR_INDICES, sprintf("SRI/SRI3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR)),
    SMA3  = file.path(DIR_INDICES, sprintf("SMA/SMA3_%d_drought_severity_sum_leq-1.0.tif", TARGET_YEAR))
  ),
  drought_leq1p5 = list(
    SPI3  = file.path(DIR_INDICES, sprintf("SPI/SPI3_%d_drought_severity_sum_leq-1.5.tif", TARGET_YEAR))
  ),


  # Output directories
  out_dir_augmented = file.path(DIR_GI_IN, "augmented"),
  out_dir_ready     = file.path(DIR_GI_IN, "ready"),
  out_dir_reports   = file.path(DIR_GI_IN, "reports"),

  # Regex patterns for identifying layer types to DROP
  name_patterns = list(
    plhiv = "^PLHIV_d[0-9]+$",
    pop   = "^POP_d[0-9]+$"
  ),

  # Regex for drought SEVERITY layers to flip (negative -> positive)
  drought_flip_pattern = "^(SPI3|SPEI3|SRI3|SMA3)_(leq1p0|leq1p5)$",

  batch_size_layers = 300
)

# Create output directories
dir.create(config$out_dir_augmented, showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_ready,     showWarnings = FALSE, recursive = TRUE)
dir.create(config$out_dir_reports,   showWarnings = FALSE, recursive = TRUE)

# Terra memory settings
dir.create(DIR_TMP, showWarnings = FALSE)
terraOptions(progress = 1)

message("Building Gi* input stacks for: ",
        paste(names(config$stacks), collapse = ", "))

# ============================================================
# 0. Utility functions
# ============================================================

vec_close <- function(a, b, tol = 1e-9) isTRUE(all(abs(a - b) <= tol))

grid_equal <- function(a, b, tol = 1e-9) {
  stopifnot(inherits(a, "SpatRaster"), inherits(b, "SpatRaster"))
  same_crs <- tryCatch(terra::same.crs(a, b), error = function(e) NA)
  if (is.na(same_crs)) same_crs <- (terra::crs(a, proj = TRUE) == terra::crs(b, proj = TRUE))
  same_res <- vec_close(terra::res(a), terra::res(b), tol)
  same_rc  <- (nrow(a) == nrow(b)) && (ncol(a) == ncol(b))
  ea <- terra::ext(a); eb <- terra::ext(b)
  same_ext <- vec_close(c(ea$xmin, ea$xmax, ea$ymin, ea$ymax),
                        c(eb$xmin, eb$xmax, eb$ymin, eb$ymax), tol)
  out <- list(crs = same_crs, res = same_res, rowcol = same_rc, extent = same_ext)
  out$all <- isTRUE(out$crs && out$res && out$rowcol && out$extent)
  out
}

report_one <- function(label, r_path, tpl, tol = 1e-9) {
  if (!file.exists(r_path)) {
    return(data.frame(label = label, path = r_path, exists = FALSE,
                      crs = NA, res = NA, rowcol = NA, extent = NA, all = NA))
  }
  r <- rast(r_path)
  ch <- grid_equal(r, tpl, tol)
  data.frame(label = label, path = r_path, exists = TRUE,
             crs = isTRUE(ch$crs), res = isTRUE(ch$res), rowcol = isTRUE(ch$rowcol),
             extent = isTRUE(ch$extent), all = isTRUE(ch$all))
}

check_alignment <- function(template_stack_path, cfg, group_name, tol = 1e-9) {
  message("\n  Alignment verification")
  message("  Template: ", basename(template_stack_path))
  tpl <- rast(template_stack_path)[[1]]

  items <- c(
    cfg$rainfall,
    cfg$heatwaves,
    setNames(cfg$drought_leq1p0, paste0(names(cfg$drought_leq1p0), "_leq1p0")),
    setNames(cfg$drought_leq1p5, paste0(names(cfg$drought_leq1p5), "_leq1p5"))
  )
  labs  <- names(items)
  paths <- unlist(items, use.names = FALSE)
  reps <- do.call(rbind, lapply(seq_along(paths), function(i) report_one(labs[i], paths[i], tpl, tol)))

  csv_file <- file.path(cfg$out_dir_reports, paste0("alignment_checks_", group_name, ".csv"))
  try({
    write.csv(reps, csv_file, row.names = FALSE)
    message("  Alignment report: ", basename(csv_file))
  }, silent = TRUE)

  n_perfect <- sum(reps$all, na.rm = TRUE)
  n_total <- sum(reps$exists, na.rm = TRUE)
  message(sprintf("  Alignment: %d/%d files perfectly aligned", n_perfect, n_total))

  misaligned <- reps[reps$exists & !reps$all, ]
  if (nrow(misaligned) > 0) {
    message("  Files needing resampling:")
    print(misaligned[, c("label", "crs", "res", "rowcol", "extent")])
  } else {
    message("  All files perfectly aligned")
  }

  reps
}

load_and_mask_layer <- function(in_path, tpl, mask_res, needs_resample = FALSE,
                                method = "bilinear", fill_zero_inside = FALSE) {
  r <- rast(in_path)

  if (needs_resample) {
    ch <- grid_equal(r, tpl)
    if (!isTRUE(ch$all)) {
      message("    Resampling: ", basename(in_path))
      r <- resample(r, tpl, method = method)
    }
  } else {
    r <- crop(r, ext(tpl))
  }

  r_m <- mask(r, mask_res)

  if (isTRUE(fill_zero_inside)) {
    v <- values(r_m)
    inside <- !is.na(values(mask_res))
    v[inside & is.na(v)] <- 0
    values(r_m) <- v
  }

  r_m
}

build_ewe_raw_layers <- function(tpl, mask_res, cfg) {
  message("  Building 9 EWE layers...")
  ras <- list()

  # -- Rainfall (2 layers, pre-aligned) --
  ras$R95pTOT <- load_and_mask_layer(cfg$rainfall$R95pTOT, tpl, mask_res, needs_resample = FALSE)
  ras$R99pTOT <- load_and_mask_layer(cfg$rainfall$R99pTOT, tpl, mask_res, needs_resample = FALSE)

  # -- Heatwaves (2 layers, pre-aligned) --
  ras$TX90p_HWMF <- load_and_mask_layer(cfg$heatwaves$TX90p_HWMF, tpl, mask_res, needs_resample = FALSE)
  ras$TX95p_HWMF <- load_and_mask_layer(cfg$heatwaves$TX95p_HWMF, tpl, mask_res, needs_resample = FALSE)

  # -- Drought severity sums leq1p0 (4 layers, pre-aligned) --
  for (nm in names(cfg$drought_leq1p0)) {
    nm_out <- paste0(nm, "_leq1p0")
    ras[[nm_out]] <- load_and_mask_layer(cfg$drought_leq1p0[[nm]], tpl, mask_res, needs_resample = FALSE)
  }

  # -- Drought severity sum leq1p5   SPI3 only (1 layer, pre-aligned) --
  for (nm in names(cfg$drought_leq1p5)) {
    nm_out <- paste0(nm, "_leq1p5")
    ras[[nm_out]] <- load_and_mask_layer(cfg$drought_leq1p5[[nm]], tpl, mask_res, needs_resample = FALSE)
  }

  n_ewe <- length(ras)
  stopifnot(n_ewe == 9)

  if (!all(vapply(ras, inherits, logical(1), "SpatRaster"))) {
    bad <- names(ras)[!vapply(ras, inherits, logical(1), "SpatRaster")]
    stop("Indicator layer did not load as a raster: ", paste(bad, collapse = ", "))
  }

  nmz <- names(ras)
  ewe_stack <- ras[[1]]
  if (length(ras) > 1) {
    for (i in 2:length(ras)) ewe_stack <- c(ewe_stack, ras[[i]])
  }
  names(ewe_stack) <- nmz

  if (nlyr(ewe_stack) != 9) {
    stop(sprintf("Expected 9 EWE layers, got %d", nlyr(ewe_stack)))
  }

  message("  9 EWE layers ready")
  ewe_stack
}

# ============================================================
# 1. Augmentation
# Drops PLHIV + POP, keeps PREV draws + country + 9 EWE
# ============================================================

augment_stack <- function(stack_path, cfg, group_name) {
  message("\n  ", strrep("-", 70))
  message("  Augmenting: ", basename(stack_path))
  message("  ", strrep("-", 70))

  if (!file.exists(stack_path)) stop("Monte Carlo draw stack not found, run step 12 first: ", stack_path)
  stk <- rast(stack_path)
  tpl <- stk[[1]]
  nms <- names(stk)

  # Drop PLHIV and POP layers
  drop_idx <- integer(0)
  if (!is.null(cfg$name_patterns$plhiv)) drop_idx <- c(drop_idx, grep(cfg$name_patterns$plhiv, nms))
  if (!is.null(cfg$name_patterns$pop))   drop_idx <- c(drop_idx, grep(cfg$name_patterns$pop,   nms))
  drop_idx <- unique(drop_idx)
  keep_idx <- setdiff(seq_along(nms), drop_idx)
  stk_keep <- stk[[keep_idx]]

  n_plhiv <- length(grep(cfg$name_patterns$plhiv, nms))
  n_pop   <- length(grep(cfg$name_patterns$pop, nms))
  base_n  <- nlyr(stk_keep)

  message(sprintf("  Removed: %d PLHIV + %d POP layers", n_plhiv, n_pop))
  message(sprintf("  Kept: %d layers (PREV draws + country)", base_n))

  # Build mask
  # Domain mask taken from the modelled PLHIV footprint;
  # 1 inside the domain, NA outside.
  .pl_idx <- grep(cfg$name_patterns$plhiv, names(stk))
  if (length(.pl_idx) == 0) stop("The draw stack has no PLHIV layers, so the analysis domain cannot be derived; rerun step 12.")
  mask_res <- terra::ifel(is.na(stk[[.pl_idx[1]]]), NA, 1L)

  # Build 9 EWE layers
  ewe_stack <- build_ewe_raw_layers(tpl, mask_res, cfg)

  # Combine (EWE only)
  out_stack <- c(stk_keep, ewe_stack)
  total <- nlyr(out_stack)
  added <- 9
  message(sprintf("  Layers: %d base + %d added = %d total", base_n, added, total))

  # Write
  out_file <- file.path(cfg$out_dir_augmented, paste0(
    tools::file_path_sans_ext(basename(stack_path)),
    "_prevOnly_EWE9.tif"
  ))

  message("  Writing augmented file...")
  writeRaster(out_stack, out_file, overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=LZW", "PREDICTOR=3", "ZLEVEL=6", "NUM_THREADS=ALL_CPUS"))
  message(sprintf("  Wrote: %s (%.2f GB)", basename(out_file), file.size(out_file) / 1024^3))

  rm(stk, stk_keep, out_stack, ewe_stack)
  gc(verbose = FALSE)

  out_file
}

# ============================================================
# 2. Transformation, written to disk in batches
#
#   sign(x) * |x|^(1/3)
#
# Gi* assumes an approximately symmetric surface. The EWE indicators are
# heavily zero-inflated, since many cells have no heatwave and no drought, and
# on that shape asinh reverses the skew instead of reducing it. The cube root
# is applied per cell, so it changes the shape of the distribution and
# therefore the z-scores.
#
# The country layer is left untransformed because it holds categorical codes.
# ============================================================

transform_to_final <- function(augmented_file, cfg, group_name) {
  message("\n  ", strrep("-", 70))
  message("  Transforming: ", basename(augmented_file))
  message("  ", strrep("-", 70))

  t_start <- Sys.time()

  # --- STEP 1: Load and identify layers ---
  r <- rast(augmented_file)
  nms <- names(r)
  n_layers <- nlyr(r)
  message(sprintf("    Total layers: %d", n_layers))

  country_idx <- which(nms == "country")
  if (length(country_idx) == 0) stop("The stack has no country layer; step 12 should have appended it")
  message(sprintf("    Country layer at index: %d (preserved untouched)", country_idx))

  to_transform_idx <- setdiff(seq_len(n_layers), country_idx)
  message(sprintf("    Layers to transform: %d", length(to_transform_idx)))

  # --- STEP 2: Identify layers to flip ---

  data_names <- nms[to_transform_idx]

  # Drought SEVERITY layers (multiply by -1: negative severity -> positive)
  drought_flip_idx <- grep(cfg$drought_flip_pattern, data_names, perl = TRUE)
  if (length(drought_flip_idx) > 0) {
    message(sprintf("    Flipping %d drought severity layers:", length(drought_flip_idx)))
    for (di in drought_flip_idx) message(sprintf("      %s", data_names[di]))
  } else {
    message("    Warning: no drought severity layers matched for negation")
  }


  # --- STEP 3: Batch processing (sign flip -> cuberoot) ---
  message("    Transform: sign(x) * |x|^(1/3)")

  batch_sz <- cfg$batch_size_layers
  n_data <- length(to_transform_idx)
  n_batches <- ceiling(n_data / batch_sz)

  batch_dir <- file.path(dirname(augmented_file),
                         paste0(".tmp_transform_", group_name))
  dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)

  batch_files <- character(n_batches)

  for (b in seq_len(n_batches)) {
    b_from <- (b - 1) * batch_sz + 1
    b_to   <- min(b * batch_sz, n_data)
    b_idx  <- b_from:b_to
    nb     <- length(b_idx)

    message(sprintf("    Batch %d/%d: layers %d-%d (%d layers)",
                    b, n_batches, b_from, b_to, nb))

    # Extract this batch of data layers
    global_idx <- to_transform_idx[b_idx]
    batch_data <- r[[global_idx]]
    batch_names <- nms[global_idx]

    # Which of THIS BATCH's local indices need flipping?
    local_drought_flip  <- which(b_idx %in% drought_flip_idx)

    # Apply drought severity sign flip
    if (length(local_drought_flip) > 0) {
      for (li in local_drought_flip) {
        batch_data[[li]] <- -1 * batch_data[[li]]
      }
      message(sprintf("      Flipped %d drought severity layers", length(local_drought_flip)))
    }


    # ============================================================
    # CUBEROOT TRANSFORMATION
    # ============================================================
    batch_data <- sign(batch_data) * abs(batch_data)^(1/3)
    names(batch_data) <- batch_names

    # Write batch to disk
    batch_file <- file.path(batch_dir, sprintf("batch_%03d.tif", b))
    writeRaster(batch_data, batch_file, overwrite = TRUE,
                datatype = "FLT4S",
                gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"))
    batch_files[b] <- batch_file

    rm(batch_data)
    gc(verbose = FALSE)
  }

  message("    Cuberoot transformation complete (all batches)")

  # --- STEP 4: Reassemble full stack and write ---

  batch_rasts <- lapply(batch_files, rast)
  data_combined <- do.call(c, batch_rasts)

  # Rebuild: insert country layer at its original position
  country_layer <- r[[country_idx]]

  parts <- vector("list", n_layers)
  k <- 1
  for (i in seq_len(n_layers)) {
    if (i == country_idx) {
      parts[[i]] <- country_layer
    } else {
      parts[[i]] <- data_combined[[k]]
      k <- k + 1
    }
  }
  stack_out <- do.call(c, parts)
  names(stack_out) <- nms

  # Write final output
  base_name <- tools::file_path_sans_ext(basename(augmented_file))
  base_name <- sub("_prevOnly_EWE9$", "", base_name)

  out_file <- file.path(cfg$out_dir_ready, paste0(base_name, "_GiStar_ready.tif"))
  message(sprintf("    Writing: %s", basename(out_file)))
  writeRaster(stack_out, out_file, overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=DEFLATE", "ZLEVEL=4", "PREDICTOR=3", "NUM_THREADS=ALL_CPUS"))

  # Clean up temp batches
  unlink(batch_dir, recursive = TRUE)

  t_end <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "mins"))

  message(sprintf("    Complete in %.1f minutes (%.2f GB)",
                  elapsed, file.size(out_file) / 1024^3))
  message("    Ready for Gi* analysis")

  rm(r, data_combined, stack_out, parts, country_layer, batch_rasts)
  gc(verbose = FALSE)

  out_file
}

# ============================================================
# Main processing loop
# ============================================================

total_files <- length(config$stacks)
completed <- 0
failed <- character(0)
outputs <- list()

message("Processing.")

for (group_name in names(config$stacks)) {
  completed <- completed + 1

  message(sprintf("# FILE %d/%d: %s", completed, total_files, group_name))

  stack_path <- config$stacks[[group_name]]

  tryCatch({

    if (!file.exists(stack_path)) {
      stop("Monte Carlo draw stack not found, run step 12 first: ", stack_path)
    }

    # Alignment check
    check_alignment(stack_path, config, group_name)

    # Augmentation
    aug_file <- augment_stack(stack_path, config, group_name)

    # Transformation
    final_file <- transform_to_final(aug_file, config, group_name)


    outputs[[group_name]] <- list(
      augmented = aug_file,
      ready     = final_file
    )

    message(sprintf("  %s complete", group_name))
    message(sprintf("    Augmented:  %s", basename(aug_file)))
    message(sprintf("    Gi*-ready:  %s", basename(final_file)))

  }, error = function(e) {
    message("\n  Error processing ", group_name, ": ", e$message)
    failed <<- c(failed, group_name)
  })

  if (completed < total_files) {
    message(sprintf("\n  PROGRESS: %d/%d complete | %d remaining",
                    completed, total_files, total_files - completed))
  }
}

message("\nAugmented stack: ", config$out_dir_augmented)
message("Gi*-ready stack: ", config$out_dir_ready)
if (length(failed)) {
  message("Failed: ", paste(failed, collapse = ", "))
} else {
  message("Done. Next: step 17.")
}
