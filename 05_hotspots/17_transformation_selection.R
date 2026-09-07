# ============================================================
# Transformation selection
# ------------------------------------------------------------
# Reports the resulting skewness of seven candidate transformations for each
# input surface, and how the one actually applied compares.
# ============================================================

source(file.path("R", "01_setup", "01_config.R"))


# Seven candidate transformations are applied to each input surface and the
# resulting skewness is reported:
#
#   raw            no transformation
#   asinh          inverse hyperbolic sine
#   log1p          log(x + c), with c set from the minimum
#   sqrt_shift     sqrt(x + c), shifted for negative values
#   rank           percentile rank on (0, 1)
#   winsor_asinh   winsorised at the 1st and 99th percentile, then asinh
#   cuberoot       sign-preserving cube root
#
# Output: the full grid of transformations by variable, the transformation
# with the smallest worst-case absolute skewness per variable, and comparison
# histograms.
#

suppressPackageStartupMessages({
  library(terra)
  library(moments)
})

# ============================================================
# CONFIG
# ============================================================

config <- list(

  # The augmented stack, before transformation.
  input_stack = file.path(DIR_GI_IN, sprintf("augmented/ALL_15_59_%d_draws_0p25_withCountry_prevOnly_EWE9.tif", TARGET_YEAR)),

  # Variables to test
  ewe_layers = c(
    "R95pTOT", "R99pTOT",
    "TX90p_HWMF", "TX95p_HWMF",
    "SPI3_leq1p0", "SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0",
    "SPI3_leq1p5"
  ),

  prev_pattern = "^PREVpct_d[0-9]+$",
  prev_n_sample = 5,  # test on 5 draws, report average

  # Layers that need sign flip BEFORE transformation
  # (drought severity: negative = worse, so multiply by -1)
  flip_layers = c(
    "SPI3_leq1p0", "SPEI3_leq1p0", "SRI3_leq1p0", "SMA3_leq1p0",
    "SPI3_leq1p5"
  ),

  # Winsorisation percentiles
  winsor_lower = 0.01,
  winsor_upper = 0.99,

  # Output
  out_dir = file.path(DIR_GI_IN, "transformation"),

  # Each transformation is assessed by the absolute skewness of the surface it
  # produces.
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
terraOptions(progress = 0)

# ============================================================
# Transformation functions
# Each takes a numeric vector (already sign-flipped if needed)
# Returns a named list: values, label, note
# ============================================================

tf_raw <- function(vals) {
  list(values = vals, label = "raw", note = "No transformation")
}

tf_asinh <- function(vals) {
  list(values = asinh(vals), label = "asinh", note = "Inverse hyperbolic sine")
}

tf_log1p <- function(vals) {
  min_val <- min(vals, na.rm = TRUE)
  if (min_val <= 0) {
    shift <- abs(min_val) + 0.001
    out <- log(vals + shift)
    note <- sprintf("log(x + %.4f)", shift)
  } else {
    out <- log1p(vals)
    note <- "log(1 + x)"
  }
  list(values = out, label = "log1p", note = note)
}

tf_sqrt_shift <- function(vals) {
  min_val <- min(vals, na.rm = TRUE)
  if (min_val < 0) {
    shift <- abs(min_val) + 0.001
    out <- sqrt(vals + shift)
    note <- sprintf("sqrt(x + %.4f)", shift)
  } else {
    out <- sqrt(vals + 0.001)
    note <- "sqrt(x + 0.001)"
  }
  list(values = out, label = "sqrt", note = note)
}

tf_rank <- function(vals) {
  out <- rank(vals, na.last = "keep", ties.method = "average")
  out <- out / max(out, na.rm = TRUE)
  list(values = out, label = "rank", note = "Percentile rank [0,1]")
}

tf_winsor_asinh <- function(vals, lower = 0.01, upper = 0.99) {
  q <- quantile(vals, probs = c(lower, upper), na.rm = TRUE)
  capped <- pmax(pmin(vals, q[2]), q[1])
  out <- asinh(capped)
  note <- sprintf("Winsorise [%.0f%%, %.0f%%] then asinh", lower * 100, upper * 100)
  list(values = out, label = "winsor_asinh", note = note)
}

tf_cuberoot <- function(vals) {
  out <- sign(vals) * abs(vals)^(1/3)
  list(values = out, label = "cuberoot", note = "Sign-preserving cube root")
}

# List of all transformations to test
all_transforms <- list(tf_raw, tf_asinh, tf_log1p, tf_sqrt_shift,
                       tf_rank, tf_winsor_asinh, tf_cuberoot)
transform_names <- c("raw", "asinh", "log1p", "sqrt", "rank", "winsor_asinh", "cuberoot")

# ============================================================
# Evaluate one variable across all transformations
# ============================================================

evaluate_variable <- function(vals_raw, layer_name, needs_flip = FALSE) {

  # Apply sign flip if needed (before any transformation)
  if (needs_flip) {
    vals_work <- -1 * vals_raw
  } else {
    vals_work <- vals_raw
  }

  results <- list()

  for (i in seq_along(all_transforms)) {
    tf_func <- all_transforms[[i]]
    tf_name <- transform_names[i]

    # Apply transformation
    tryCatch({
      if (tf_name == "winsor_asinh") {
        tf_out <- tf_func(vals_work, config$winsor_lower, config$winsor_upper)
      } else {
        tf_out <- tf_func(vals_work)
      }

      tv <- tf_out$values
      tv <- tv[is.finite(tv)]

      if (length(tv) < 10) {
        results[[i]] <- data.frame(
          layer = layer_name, transform = tf_name,
          note = tf_out$note,
          skewness = NA, abs_skewness = NA,
          kurtosis = NA, mean = NA, sd = NA,
          min = NA, max = NA, n_valid = length(tv),
          sign_flipped = needs_flip,
          stringsAsFactors = FALSE
        )
        next
      }

      skew <- skewness(tv)
      kurt <- kurtosis(tv)

      results[[i]] <- data.frame(
        layer         = layer_name,
        transform     = tf_name,
        note          = tf_out$note,
        skewness      = round(skew, 4),
        abs_skewness  = round(abs(skew), 4),
        kurtosis      = round(kurt, 4),
        mean          = round(mean(tv), 4),
        sd            = round(sd(tv), 4),
        min           = round(min(tv), 4),
        max           = round(max(tv), 4),
        n_valid       = length(tv),
        sign_flipped  = needs_flip,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      results[[i]] <<- data.frame(
        layer = layer_name, transform = tf_name,
        note = paste("ERROR:", e$message),
        skewness = NA, abs_skewness = NA,
        kurtosis = NA, mean = NA, sd = NA,
        min = NA, max = NA, n_valid = 0,
        sign_flipped = needs_flip,
        stringsAsFactors = FALSE
      )
    })
  }

  do.call(rbind, results)
}

# ============================================================
# MAIN EXECUTION
# ============================================================

cat("\n")
cat("  Transformation selection grid\n")
cat(strrep("=", 80), "\n\n")

# Load stack
cat(sprintf("  Loading: %s\n", basename(config$input_stack)))
r <- rast(config$input_stack)
all_names <- names(r)
cat(sprintf("  Layers: %d\n\n", nlyr(r)))

all_grid_results <- list()

# -- Process EWE layers --
cat("  Indicator layers\n")

for (lyr in config$ewe_layers) {
  if (!lyr %in% all_names) {
    cat(sprintf("    %-25s not found, skipped\n", lyr))
    next
  }

  needs_flip <- lyr %in% config$flip_layers
  vals <- values(r[[lyr]], mat = FALSE)
  vals <- vals[is.finite(vals)]

  result <- evaluate_variable(vals, lyr, needs_flip)
  all_grid_results[[lyr]] <- result

  # Find best
  best_idx <- which.min(result$abs_skewness)
  best <- result[best_idx, ]

  cat(sprintf("    %-25s best = %-15s |skew| = %.3f  (raw = %.3f)%s\n",
              lyr, best$transform, best$abs_skewness,
              result$abs_skewness[result$transform == "raw"],
              ifelse(needs_flip, "  [flipped]", "")))

}

# -- Process PREV draws (sample, then average) --
cat("\n  HIV PREVALENCE DRAWS\n")

prev_idx <- grep(config$prev_pattern, all_names)
cat(sprintf("    Total PREV draws: %d\n", length(prev_idx)))

if (length(prev_idx) > 0) {
  if (length(prev_idx) > config$prev_n_sample) {
    sample_idx <- prev_idx[round(seq(1, length(prev_idx),
                                      length.out = config$prev_n_sample))]
  } else {
    sample_idx <- prev_idx
  }

  # Collect skewness for each transform across sampled draws
  prev_skew_matrix <- matrix(NA, nrow = length(sample_idx),
                              ncol = length(transform_names))
  colnames(prev_skew_matrix) <- transform_names

  for (si in seq_along(sample_idx)) {
    idx <- sample_idx[si]
    lyr_name <- all_names[idx]
    vals <- values(r[[idx]], mat = FALSE)
    vals <- vals[is.finite(vals)]

    result <- evaluate_variable(vals, lyr_name, needs_flip = FALSE)
    prev_skew_matrix[si, ] <- result$abs_skewness
  }

  # Average across draws
  avg_abs_skew <- colMeans(prev_skew_matrix, na.rm = TRUE)
  best_tf <- transform_names[which.min(avg_abs_skew)]

  cat(sprintf("    Sampled %d draws, averaged |skewness| per transform:\n",
              length(sample_idx)))
  for (i in seq_along(transform_names)) {
    marker <- ifelse(transform_names[i] == best_tf, " <-- BEST", "")
    cat(sprintf("      %-15s avg |skew| = %.3f%s\n",
                transform_names[i], avg_abs_skew[i], marker))
  }

  # Store as summary result
  prev_summary <- data.frame(
    layer = "PREV_draws_average",
    transform = transform_names,
    note = "Average across sampled draws",
    skewness = NA,   
    abs_skewness = round(avg_abs_skew, 4),
    kurtosis = NA,
    mean = NA, sd = NA, min = NA, max = NA,
    n_valid = length(sample_idx),
    sign_flipped = FALSE,
    stringsAsFactors = FALSE
  )
  all_grid_results[["PREV_draws_average"]] <- prev_summary

}

rm(r); gc(verbose = FALSE)

# ============================================================
# Compile full grid and recommendations
# ============================================================

full_grid <- do.call(rbind, all_grid_results)
rownames(full_grid) <- NULL

# Best per variable
best_per_var <- do.call(rbind, lapply(split(full_grid, full_grid$layer), function(df) {
  best_idx <- which.min(df$abs_skewness)
  df[best_idx, ]
}))
best_per_var <- best_per_var[order(best_per_var$abs_skewness), ]
rownames(best_per_var) <- NULL

APPLIED <- "cuberoot"
applied_per_var <- do.call(rbind, lapply(split(full_grid, full_grid$layer), function(df) {
  hit <- df[df$transform == APPLIED, ]
  if (nrow(hit) >= 1) hit[1, ] else df[which.min(df$abs_skewness), ]
}))
rownames(applied_per_var) <- NULL

runner_up_per_var <- do.call(rbind, lapply(split(full_grid, full_grid$layer), function(df) {
  df_sorted <- df[order(df$abs_skewness), ]
  if (nrow(df_sorted) >= 2) df_sorted[2, ] else df_sorted[1, ]
}))
rownames(runner_up_per_var) <- NULL

# ============================================================
# Display results
# ============================================================

cat("\n\n")
cat("  Full transformation grid\n")
cat(strrep("=", 80), "\n\n")

# Print pivot-style table
layers_order <- c(config$ewe_layers, "PREV_draws_average")

cat(sprintf("  %-25s", "Layer"))
for (tn in transform_names) cat(sprintf("  %12s", tn))
cat("   best\n")
cat(sprintf("  %s\n", strrep("-", 25 + length(transform_names) * 14 + 20)))

for (lyr in layers_order) {
  lyr_data <- full_grid[full_grid$layer == lyr, ]
  if (nrow(lyr_data) == 0) next

  cat(sprintf("  %-25s", lyr))

  min_skew <- min(lyr_data$abs_skewness, na.rm = TRUE)

  for (tn in transform_names) {
    row <- lyr_data[lyr_data$transform == tn, ]
    if (nrow(row) > 0 && !is.na(row$abs_skewness)) {
      marker <- ifelse(row$abs_skewness == min_skew, "*", " ")
      cat(sprintf("  %11.3f%s", row$abs_skewness, marker))
    } else {
      cat(sprintf("  %12s", "N/A"))
    }
  }

  best_row <- lyr_data[which.min(lyr_data$abs_skewness), ]
  cat(sprintf("   %s", best_row$transform))
  cat("\n")
}

cat("\n  (* = best for that variable)\n")

# ============================================================
# Recommendation table
# ============================================================

cat("\n\n")
cat("  Candidate transformations per variable\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("  %-25s  %-13s %8s  %-13s %8s  %8s\n",
            "Layer", "Best", "|Skew|", "Runner-up", "|Skew|", APPLIED))
cat(sprintf("  %s\n", strrep("-", 84)))

for (lyr in layers_order) {
  b  <- best_per_var[best_per_var$layer == lyr, ]
  ru <- runner_up_per_var[runner_up_per_var$layer == lyr, ]
  ap <- applied_per_var[applied_per_var$layer == lyr, ]

  if (nrow(b) > 0) {
    cat(sprintf("  %-25s  %-13s %8.3f  %-13s %8.3f  %8.3f\n",
                lyr, b$transform, b$abs_skewness,
                ifelse(nrow(ru) > 0, ru$transform, "N/A"),
                ifelse(nrow(ru) > 0, ru$abs_skewness, NA),
                ifelse(nrow(ap) > 0, ap$abs_skewness, NA)))
  }
}


# ============================================================
# SAVE OUTPUTS
# ============================================================

# Full grid
csv_grid <- file.path(config$out_dir, "transformation_grid_full.csv")
write.csv(full_grid, csv_grid, row.names = FALSE)

# Recommendations
csv_rec <- file.path(config$out_dir, "transformation_recommendations.csv")
rec_out <- merge(
  best_per_var[, c("layer", "transform", "abs_skewness", "sign_flipped")],
  runner_up_per_var[, c("layer", "transform", "abs_skewness")],
  by = "layer", suffixes = c("_best", "_runnerup")
)

rec_out <- merge(rec_out,
                 setNames(applied_per_var[, c("layer", "abs_skewness")],
                          c("layer", paste0("abs_skewness_", APPLIED))),
                 by = "layer")
rec_out <- rec_out[order(rec_out$layer), ]
write.csv(rec_out, csv_rec, row.names = FALSE)

message("Done. Next: step 18.")
