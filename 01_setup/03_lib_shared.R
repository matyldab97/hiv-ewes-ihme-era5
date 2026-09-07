# ============================================================
# Helpers shared by the exposure and hotspot steps
# ------------------------------------------------------------
# Country assignment, classification breaks, country-layer lookups, number
# formatting, the two ggplot themes, and the PLHIV counters that read the draw matrix.
#
# The counters expect two objects in the calling script: `plhiv_mat`, the
# cells-by-draws matrix the intervals are computed from, and `plhiv_med`, one
# value per cell, which is what decides whether a cell carries population at
# all. get_label() expects `config$var_labels`.
#
# Sourced by steps 13 to 19. Not run directly.
# ============================================================

fmt3 <- function(x) sprintf("%.3f", x)


fmt_num <- function(x) {
  ifelse(x >= 1e6, sprintf("%.1fM", x / 1e6),
  ifelse(x >= 1e3, sprintf("%.0fK", x / 1e3),
                   sprintf("%.0f",  x)))
}

is_discrete <- function(v) {
  length(v) > 0 && all(v == round(v)) && (max(v) - min(v)) < 50
}

classify_vals <- function(vals, br) {
  if (is.null(br) || length(br) < 2) return(rep(NA_integer_, length(vals)))
  as.integer(cut(vals, breaks = br, include.lowest = TRUE, right = FALSE, labels = FALSE))
}

compute_iqr_breaks <- function(vp) {
  if (length(vp) == 0) return(NULL)
  br <- as.numeric(quantile(vp, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  if (is_discrete(vp)) br <- round(br)
  sort(unique(br))
}

compute_range_breaks <- function(vp) {
  if (length(vp) == 0) return(NULL)
  vmin <- min(vp); vmax <- max(vp)
  if (vmin == vmax) return(c(vmin - 0.5, vmax + 0.5))
  discrete <- is_discrete(vp)
  if (discrete) {
    imin <- as.integer(vmin); imax <- as.integer(vmax)
    if ((imax - imin + 1) <= 4) br <- c(imin - 0.5, seq(imin, imax) + 0.5)
    else br <- seq(imin - 0.5, imax + 0.5, length.out = 5)
  } else { br <- seq(vmin, vmax, length.out = 5) }
  sort(unique(br))
}

# A country map is drawn by joining results to the boundary polygons by name.
# The join comes back empty when the draw stack carries country codes without
# their labels, which happens when a GeoTIFF is moved without its .aux.xml
# sidecar. 
drawable_or_skip <- function(map_sf, panel_col, what) {
  d <- map_sf[!is.na(map_sf[[panel_col]]), ]
  if (nrow(d) == 0) {
    message("  Skipped: ", what,
            " (no country matched between the draw stack and the boundary file)")
    return(NULL)
  }
  d
}

find_country_layer <- function(stk) {
  idx <- which(tolower(names(stk)) == "country")
  if (length(idx) == 1) return(stk[[idx]])
  for (i in seq_len(nlyr(stk))) {
    lev <- tryCatch(levels(stk[[i]]), error = function(e) NULL)
    if (is.list(lev) && length(lev) >= 1 && is.data.frame(lev[[1]]) &&
        ("country" %in% names(lev[[1]]))) return(stk[[i]])
  }
  stop("The draw stack has no country layer; step 12 appends it when it writes the stack.")
}

# The value-to-name table for the country layer. terra carries it on the raster
# itself; a GeoTIFF keeps it in the .tif.aux.xml sidecar beside the .tif, so a
# .tif copied on its own arrives without country names.
# When the layer has lost the table, the same table is read from the source
# stack instead. 
get_country_levels <- function(r_country, r_country_ref = NULL) {
  as_table <- function(r) {
    if (is.null(r)) return(NULL)
    lev <- tryCatch(levels(r), error = function(e) NULL)
    if (!is.list(lev) || length(lev) < 1 || !is.data.frame(lev[[1]]) ||
        nrow(lev[[1]]) == 0) return(NULL)
    df <- lev[[1]]
    if (!all(c("value", "country") %in% names(df))) {
      df <- df[, c(1, ncol(df)), drop = FALSE]; names(df) <- c("value", "country")
    }
    df[, c("value", "country"), drop = FALSE]
  }
  df <- as_table(r_country)
  if (is.null(df)) df <- as_table(r_country_ref)
  if (is.null(df))
    stop("The country layer carries no value-to-name table, so countries could ",
         "only be labelled by number. A GeoTIFF keeps that table in its ",
         ".tif.aux.xml sidecar: copy that file next to the .tif and run again.")
  df
}

# Name for one country code, read from the table get_country_levels() returned.
country_name <- function(cid, levels_df) {
  nm <- levels_df$country[match(cid, levels_df$value)]
  if (is.na(nm))
    stop("Country code ", cid, " has no name in the country level table. ",
         "The country layer and the table it came from are out of step; ",
         "recopy the layer and its .tif.aux.xml sidecar from step 12.")
  as.character(nm)
}

compute_country_raster <- function(tpl_rast, gadm_path, name_col = "COUNTRY") {
  countries <- st_read(gadm_path, quiet = TRUE)
  tpl_id <- rast(tpl_rast); tpl_id[] <- 1:ncell(tpl_id)
  tpl_crs_wkt <- crs(tpl_id)
  if (!identical(st_crs(countries)$wkt, tpl_crs_wkt))
    countries <- st_transform(countries, tpl_crs_wkt)
  overlaps <- exact_extract(tpl_id, countries, include_cell = TRUE, progress = TRUE)
  maj <- rep(NA_character_, ncell(tpl_id)); maxcov <- rep(0, ncell(tpl_id))
  for (i in seq_along(overlaps)) {
    df <- overlaps[[i]]
    if (!is.null(df) && nrow(df) > 0) {
      cn <- countries[[name_col]][i]
      cov <- tapply(df$coverage_fraction, df$cell, sum, na.rm = TRUE)
      for (cid in names(cov)) {
        ci <- as.integer(cid)
        if (cov[cid] > maxcov[ci]) { maj[ci] <- cn; maxcov[ci] <- cov[cid] }
      }
    }
  }
  list(assignment = maj, countries_sf = countries)
}

get_label <- function(var) {
  if (var %in% names(config$var_labels)) config$var_labels[var] else var
}

count_plhiv_ci <- function(mask) {
  if (is.null(plhiv_mat)) {
    val <- sum(plhiv_med[mask], na.rm = TRUE)
    return(data.frame(plhiv_mean = round(val), plhiv_lower = NA_real_, plhiv_upper = NA_real_))
  }
  draw_sums <- colSums(plhiv_mat[mask, , drop = FALSE], na.rm = TRUE)
  data.frame(plhiv_mean = round(mean(draw_sums)),
             plhiv_lower = round(quantile(draw_sums, 0.025)),
             plhiv_upper = round(quantile(draw_sums, 0.975)))
}

pct_plhiv_ci <- function(subset_mask, total_mask) {
  if (is.null(plhiv_mat)) {
    num <- sum(plhiv_med[subset_mask], na.rm = TRUE)
    den <- sum(plhiv_med[total_mask], na.rm = TRUE)
    pct <- if (den > 0) round(100 * num / den, 2) else NA
    return(data.frame(pct_mean = pct, pct_lower = NA_real_, pct_upper = NA_real_))
  }
  num_draws <- colSums(plhiv_mat[subset_mask, , drop = FALSE], na.rm = TRUE)
  den_draws <- colSums(plhiv_mat[total_mask, , drop = FALSE], na.rm = TRUE)
  pct_draws <- ifelse(den_draws > 0, 100 * num_draws / den_draws, NA_real_)
  data.frame(pct_mean = round(mean(pct_draws, na.rm = TRUE), 2),
             pct_lower = round(quantile(pct_draws, 0.025, na.rm = TRUE), 2),
             pct_upper = round(quantile(pct_draws, 0.975, na.rm = TRUE), 2))
}

theme_map <- function(base_size = 11) {
  theme_void(base_size = base_size) +
    theme(plot.title=element_text(face="bold",size=base_size+3,hjust=0),
          plot.subtitle=element_text(colour="grey40",size=base_size,hjust=0),
          strip.text=element_text(face="bold",size=base_size+1,margin=margin(t=3,b=4)),
          strip.clip="off",
          legend.position="bottom", legend.key.width=unit(1.5,"cm"),
          plot.margin=margin(5,5,5,5))
}

theme_ewe <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(plot.title=element_text(face="bold",size=base_size+4,hjust=0),
          plot.subtitle=element_text(colour="grey40",size=base_size,hjust=0),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.y=element_line(colour="grey85", linewidth=0.3),
          panel.grid.minor.y=element_line(colour="grey92", linewidth=0.2),
          axis.title.y=element_text(size=base_size), legend.position="bottom",
          strip.text=element_text(face="bold",size=base_size+1,margin=margin(t=3,b=4)),
          strip.clip="off",
          plot.margin=margin(10,15,10,10))
}
