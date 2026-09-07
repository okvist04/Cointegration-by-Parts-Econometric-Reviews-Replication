# ============================================================
# run_scan_test.R
#
# Apply FIEG(r0), BIEG(r0,1), or GIEG(r0) to a real dataset, and
# (optionally) attach significance against the critical values fitted
# in fit_response_surface.R.
#
# Requires generate_cv.R to have been sourced (for scan_full_cpp(), via
# adf_scan_fast.cpp) and, if you want critical values attached,
# fit_response_surface.R too. plot.scan_test() requires ggplot2.
# ============================================================

source("R/generate_cv.R")
source("R/fit_response_surface.R")
library(ggplot2)

.scheme_to_code <- c(FIEG = "fw", BIEG = "bw", GIEG = "all")

# ------------------------------------------------------------
# run_scan_test
# ------------------------------------------------------------
#
# y, X       : the real data. X is coerced to a matrix (n_regressors
#              columns); y and X must have the same number of rows.
# scheme     : "FIEG", "BIEG", or "GIEG" -- matches the paper's notation.
# r0         : the trimming fraction. MUST match whatever r0 the critical
#              values / response surface you'll compare against were
#              calibrated with -- there is no way to check this from the
#              data alone, so get this right by construction (e.g. define
#              it once as a shared constant alongside run_generation_cv.R's
#              configuration, rather than retyping the number here).
# include_intercept, max_lag : must also match the calibration settings
#              (max_lag = NULL reproduces the same automatic
#              12*(n/100)^(1/4) rule used during calibration).
# return_full_scan : if TRUE (default), also returns the full window-by-
#              window scan path in $scan, e.g. for plotting EG(r) (FIEG,
#              BIEG) or exporting the EG(r1,r2) surface (GIEG).
#
# Returns an object of class "scan_test" with the statistic and the
# location of the extremizing window, in both observation-index and
# sample-fraction (r) terms.

run_scan_test <- function(
    y,
    X,
    scheme = c("FIEG", "BIEG", "GIEG"),
    r0,
    include_intercept = TRUE,
    max_lag = NULL,
    return_full_scan = TRUE,
    dates = NULL,
    min_window = NULL
) {
  
  scheme <- match.arg(scheme)
  scan_type <- .scheme_to_code[[scheme]]
  
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  Tt <- length(y)
  
  if (nrow(X) != Tt)
    stop("y and X must have the same number of observations.")
  if (r0 <= 0 || r0 >= 1)
    stop("r0 must be in (0,1).")
  if (!is.null(dates) && length(dates) != Tt)
    stop("dates must have length equal to the number of observations.")
  
  m_T <- floor(r0 * Tt)

  # Enforce an absolute floor independent of r0*T, since on small
  # subsamples the fraction-based window can shrink to a size where
  # the ADF regression is nearly saturated (too few residual degrees
  # of freedom relative to lags), producing spuriously extreme
  # t-statistics rather than genuine evidence of cointegration.
  if (!is.null(min_window)) {
    m_T <- max(m_T, min_window)
  }
  
  if (m_T < 10)
    stop(sprintf(
      "r0 * T = %.1f gives a minimum window of only %d observations -- too small for the ADF regression (needs >= 10). Use a larger r0 or a longer series.",
      r0 * Tt, m_T
    ))

   if (m_T >= Tt)
    stop(sprintf(
      "Minimum window (%d obs) is >= the sample size (%d obs) -- no valid window exists. Reduce min_window or drop this subsample.",
      m_T, Tt
    ))  
  
  ml <- if (is.null(max_lag)) -1L else as.integer(max_lag)
  
  res <- scan_full_cpp(
    y_r = as.numeric(y),
    X_r = X,
    scan_type = scan_type,
    min_obs = m_T,
    min_start = 1L,
    min_window_obs = m_T,
    include_intercept = include_intercept,
    max_lag = ml,
    return_full_scan = return_full_scan
  )
  
  if (is.na(res$min_t)) {
    warning("No window produced a finite statistic -- check the data (e.g. constant series, insufficient variation, or all candidate windows rank-deficient).")
  }
  
  out <- list(
    scheme = scheme,
    statistic = res$min_t,
    T = Tt,
    n_regressors = ncol(X),
    r0 = r0,
    best_start = res$best_start,
    best_end = res$best_end,
    best_lag = res$best_lag,
    r1_hat = if (is.na(res$best_start)) NA_real_ else (res$best_start - 1) / Tt,
    r2_hat = if (is.na(res$best_end)) NA_real_ else res$best_end / Tt
  )
  
  if (!is.null(dates) && !is.na(res$best_start)) {
    out$date_start <- dates[res$best_start]
    out$date_end <- dates[res$best_end]
  }
  
  if (return_full_scan) {
    out$scan <- data.table(
      start = res$start,
      end = res$end,
      tstat = res$tstat,
      lag = res$lag
    )
  }
  
  class(out) <- "scan_test"
  out
}

# ------------------------------------------------------------
# add_critical_values
# ------------------------------------------------------------
#
# Attaches critical values and reject/fail-to-reject verdicts from a
# fitted response-surface coefficient table (see fit_response_surface.R)
# to a scan_test result. Kept separate from run_scan_test() so you can
# run the scan itself before the full calibration grid is finished, and
# attach significance once coef_table is available.

add_critical_values <- function(
    test_result,
    coef_table,
    probs = c(0.01, 0.05, 0.10)
) {
  
  if (!inherits(test_result, "scan_test"))
    stop("test_result must come from run_scan_test().")
  
  scan_type <- .scheme_to_code[[test_result$scheme]]
  
  cvs <- vapply(probs, function(p) {
    get_critical_value(
      T = test_result$T,
      n_regressors = test_result$n_regressors,
      scan_type = scan_type,
      prob = p,
      coef_table = coef_table
    )
  }, numeric(1))
  
  names(cvs) <- paste0(format(100 * probs, trim = TRUE), "%")
  
  test_result$critical_values <- cvs
  # Rejection at level p: the observed statistic is more extreme (more
  # negative) than the p-quantile of the null distribution.
  test_result$reject <- setNames(test_result$statistic < cvs, names(cvs))
  
  test_result
}

# ------------------------------------------------------------
# test_cointegration
# ------------------------------------------------------------
#
# Convenience one-shot wrapper: runs the scan and, if coef_table is
# supplied, attaches critical values in the same call.

test_cointegration <- function(
    y,
    X,
    scheme = c("FIEG", "BIEG", "GIEG"),
    r0,
    coef_table = NULL,
    probs = c(0.01, 0.05, 0.10),
    include_intercept = TRUE,
    max_lag = NULL,
    return_full_scan = TRUE,
    dates = NULL,
    min_window = NULL
) {
  
  res <- run_scan_test(
    y = y, X = X, scheme = scheme, r0 = r0,
    include_intercept = include_intercept, max_lag = max_lag,
    return_full_scan = return_full_scan, dates = dates, min_window = min_window
  )
  
  if (!is.null(coef_table)) {
    res <- add_critical_values(res, coef_table, probs = probs)
  }
  
  res
}

# ------------------------------------------------------------
# print / plot methods
# ------------------------------------------------------------

print.scan_test <- function(x, ...) {
  
  cat(sprintf(
    "%s(r0 = %.3f) test statistic: %.4f  (T = %d, k = %d)\n",
    x$scheme, x$r0, x$statistic, x$T, x$n_regressors
  ))
  
  if (!is.na(x$best_start)) {
    loc <- sprintf(
      "Extremum at observations %d:%d (r1 = %.3f, r2 = %.3f), lag = %d",
      x$best_start, x$best_end, x$r1_hat, x$r2_hat, x$best_lag
    )
    if (!is.null(x$date_start)) {
      loc <- paste0(loc, sprintf(" [%s to %s]", x$date_start, x$date_end))
    }
    cat(loc, "\n")
  }
  
  if (!is.null(x$critical_values)) {
    cat("\nCritical values and rejection of H0 (no cointegration):\n")
    for (nm in names(x$critical_values)) {
      cat(sprintf(
        "  %-4s critical value: %8.4f   reject H0: %s\n",
        nm, x$critical_values[[nm]], x$reject[[nm]]
      ))
    }
  } else {
    cat("\n(No critical values attached -- pass coef_table to test_cointegration() or add_critical_values() to get a decision.)\n")
  }
  
  invisible(x)
}

# Diagnostic plot of the scan path. FIEG/BIEG scan a single endpoint, so
# EG(r) is plotted directly against r with the extremum and (if present)
# critical value marked. GIEG scans both endpoints (a 2D surface), so
# rather than guess at a visualization, this points you to $scan for a
# custom heatmap/contour instead of producing one that may not fit how
# you want to present it in the paper.

plot.scan_test <- function(x, level = "5%", ...) {
  
  if (is.null(x$scan))
    stop("No scan path stored -- re-run with return_full_scan = TRUE.")
  
  if (x$scheme == "GIEG") {
    message(
      "GIEG scans a 2D (r1, r2) surface -- no default plot is provided. ",
      "Use x$scan (columns: start, end, tstat, lag) to build a contour ",
      "or heatmap suited to how you want to present it."
    )
    return(invisible(x))
  }
  
  r <- if (x$scheme == "FIEG") x$scan$end / x$T else x$scan$start / x$T
  df <- data.frame(r = r, tstat = x$scan$tstat)
  df <- df[order(df$r), ]
  
  r_hat <- if (x$scheme == "FIEG") x$r2_hat else x$r1_hat
  has_cv <- !is.null(x$critical_values) && level %in% names(x$critical_values)
  
  # Colours are mapped through aes() (rather than passed as fixed
  # geom arguments) purely so each line gets its own legend entry --
  # the values plotted don't depend on this, it's a legend-building trick.
  colour_values <- c("EG(r)" = "black", "extremum" = "grey50")
  
  p <- ggplot(df, aes(x = .data$r, y = .data$tstat)) +
    geom_line(aes(colour = "EG(r)")) +
    geom_vline(aes(xintercept = r_hat, colour = "extremum"), linetype = "dashed")
  
  if (has_cv) {
    cv <- x$critical_values[[level]]
    cv_label <- paste0(level, " critical value")
    colour_values[cv_label] <- "firebrick"
    p <- p + geom_hline(aes(yintercept = cv, colour = cv_label), linetype = "dotted")
  }
  
  p <- p +
    scale_colour_manual(name = NULL, values = colour_values) +
    labs(
      x = "r", y = "EG(r)",
      title = sprintf("%s scan path", x$scheme),
      subtitle = sprintf("T = %d, k = %d", x$T, x$n_regressors)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  invisible(p)
}

# ============================================================
# Example usage (uncomment once you have real data and, if you want
# significance attached, a fitted coef_table from fit_response_surface.R):
#
# y <- your_y_vector
# X <- your_X_matrix_or_vector   # n_regressors columns
#
# coef_table <- fread("output/cv_output/response_surface_coefficients.csv")
#
# result <- test_cointegration(
#   y = y, X = X, scheme = "GIEG", r0 = 0.15,   # r0 MUST match calibration
#   coef_table = coef_table
# )
# print(result)
# plot(result)   # FIEG/BIEG only; GIEG points you to result$scan instead
# ============================================================