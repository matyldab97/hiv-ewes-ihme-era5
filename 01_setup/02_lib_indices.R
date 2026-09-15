# ============================================================
# Shared engine for the standardised drought indices
# ------------------------------------------------------------
# Provides the monthly reference builder, the parametric index fitter (gamma for
# SPI and SRI, log-logistic for SPEI), the empirical per-calendar-month
# transformation used by SMA, and the annual severity sum.
#
# Sourced by steps 07 to 10. Not run directly.
# ============================================================

source(file.path("01_setup", "00_paths.R"))

suppressPackageStartupMessages({ library(SPEI); library(zoo) })

dpm_of <- function(year) if (is_leap_year(year))
  c(31,29,31,30,31,30,31,31,30,31,30,31) else c(31,28,31,30,31,30,31,31,30,31,30,31)

# ---- daily -> 12 monthly layers, fun = "sum" or "mean" ------------------
# A complete year is assumed, which 00_inputs/01_era5_inputs.R verifies by
# counting the daily layers of every input file.
monthly_from_daily <- function(daily_r, year, fun = c("sum","mean")) {
  fun <- match.arg(fun)
  idx <- rep(1:12, dpm_of(year))
  tapp(daily_r, idx, if (fun == "sum") sum else mean, na.rm = TRUE)
}

# ---- monthly stack over ANALYSIS_YEARS; layers named YYYY_MM ------------------
# The default span covers the reference period and the target year. The fitting
# period stays REF_YEARS: standardised_stack() passes it as ref.start/ref.end,
# and sma_target_year() selects it by name.
build_monthly_ref <- function(var, fun, years = ANALYSIS_YEARS) {
  lst <- lapply(years, function(y) {
    mm <- monthly_from_daily(load_year(var, y), y, fun)
    names(mm) <- sprintf("%d_%02d", y, 1:12); mm
  })
  rast(lst)
}

# ---- clamp, applied to every standardised index -------------------------
# A gamma or log-logistic fit on a baseline that is mostly zeros is not
# identifiable, and the resulting index diverges rather than saturating. The
# bound keeps those cells finite so that downstream severity sums stay
# comparable.
#
# NaN is mapped to NA rather than to the bound. NaN means the fit returned no
# value at all, which is a different thing from a value that ran off the scale.
clamp_uniform <- function(x, b = CLAMP_BOUND) {
  x[is.nan(x)] <- NA_real_
  ii <- is.infinite(x); x[ii] <- sign(x[ii]) * b
  x[!is.na(x) & x >  b] <-  b
  x[!is.na(x) & x < -b] <- -b
  x
}

# ---- gamma / log-logistic standardised index ----------------------------
# kind "spi" -> SPEI::spi ; "spei" -> SPEI::spei (the stack is then P - PET)
#
#
# Two failure modes are counted and reported:
#
#   skipped   the series is all missing, constant, or all zero
#   clamped   the fit ran but the index left [-CLAMP_BOUND, CLAMP_BOUND]. On a
#             zero-inflated baseline the gamma fit is not identifiable and the
#             index diverges. This is the arid and dry-season failure that makes
#             SPI-3 unreliable over the Sahel in December to April.
#
# A clamped cell can be a fit failure carried forward as an extreme value.
standardised_stack <- function(monthly_stack, scale, kind = c("spi","spei"),
                               clamp = clamp_uniform) {
  kind    <- match.arg(kind)
  fit_fun <- if (kind == "spi") SPEI::spi else SPEI::spei

  y0  <- as.integer(substr(names(monthly_stack)[1], 1, 4))
  ref <- list(start = c(REF_YEARS[1], 1), end = c(REF_YEARS[length(REF_YEARS)], 12))
  M   <- t(as.matrix(monthly_stack))
  n   <- ncol(M)

  fit_one <- function(s) {
    if (all(is.na(s)) || sd(s, na.rm = TRUE) == 0 ||
        sum(abs(s), na.rm = TRUE) == 0)
      return(list(v = rep(NA_real_, length(s)), skipped = 1L, clamped = 0L, failed = 0L))
    raw <- tryCatch(as.numeric(fit_fun(ts(s, start = c(y0, 1), frequency = 12),
                                       scale = scale,
                                       ref.start = ref$start,
                                       ref.end = ref$end)$fitted),
                    error = function(e) NULL)
    if (is.null(raw))
      return(list(v = rep(NA_real_, length(s)), skipped = 0L, clamped = 0L, failed = 1L))
    n_clamped <- sum(is.infinite(raw) | (!is.na(raw) & abs(raw) > CLAMP_BOUND))
    list(v = clamp(raw), skipped = 0L, clamped = n_clamped, failed = 0L)
  }

  out   <- matrix(NA_real_, nrow(M), n)
  tally <- c(skipped = 0L, clamped = 0L, failed = 0L)
  for (k in seq_len(n)) {
    one <- fit_one(M[, k])
    out[, k] <- one$v
    tally <- tally + c(one$skipped, one$clamped, one$failed)
  }

  n_fitted <- n - tally[["skipped"]] - tally[["failed"]]
  message(sprintf("    %s-%d: %d cells fitted, %d skipped (flat or empty baseline), %d fit errors",
                  toupper(kind), scale, n_fitted, tally[["skipped"]], tally[["failed"]]))
  if (tally[["clamped"]] > 0) {
    total_months <- n_fitted * nrow(M)
    message(sprintf("    %s-%d: %d of %d cell-months clamped at +/-%g (%.2f%%), which marks cells where the distribution could not be fitted",
                    toupper(kind), scale, tally[["clamped"]], total_months,
                    CLAMP_BOUND, 100 * tally[["clamped"]] / max(total_months, 1)))
  }

  r <- setValues(monthly_stack, t(out)); names(r) <- names(monthly_stack); r
}

# ---- pull target-year layers from a YYYY_MM stack ------------------------------------------------------------
target_year_layers <- function(stk, year = TARGET_YEAR) {
  stk[[which(names(stk) %in% sprintf("%d_%02d", year, 1:12))]]
}

# ---- drought severity sum -----------------------------------------------
# Sums the index over the months at or below the threshold, so both duration
# and intensity contribute; a cell with no such month scores 0,
# which keeps it in the denominator as unexposed.
#
# The sum is negative by construction, because the index is negative in
# drought. Steps 13 and 16 negate it so that larger always means worse.
severity_sum <- function(year_stack, thr) {
  dv <- year_stack; dv[dv > thr] <- NA
  sev <- app(dv, sum, na.rm = TRUE)
  cnt <- app(year_stack <= thr, sum, na.rm = TRUE)
  sev[cnt == 0] <- 0
  sev
}
