# =====================================================================
# check_tail_fit.R  --  does the calibration grid reach its limit?
#
#   Rscript check_tail_fit.R
#
# The response surface's beta_inf is only interpretable as an asymptotic
# critical value if the tail of the grid is actually converging toward
# it. Two checks, using every sample size (weighted by its own Monte
# Carlo standard error), on criteria fixed in advance:
#
#   A. TAIL FIT. Standardized residuals of the largest sample sizes
#      under the full fit should have |resid| < 2.
#   B. EXTRAPOLATION GAP. beta_inf should be within 2 combined standard
#      errors of the largest simulated critical value.
#
# The grid "reaches the limit" if a MAJORITY of (case, N, test, level)
# combinations pass BOTH criteria.
#
# Environment: TABLES [results/tables]   TAIL_T [700]
# =====================================================================

source("R/response_surface.R")

TABLES <- Sys.getenv("TABLES", "results/tables")
TAIL_T <- as.integer(Sys.getenv("TAIL_T", "700"))
LEVELS <- c(0.01, 0.05, 0.10)

main <- function() {
  f <- file.path(TABLES, "quantiles.csv")
  if (!file.exists(f)) stop(sprintf("%s not found -- run make_tables.R first", f))
  q <- read.csv(f, stringsAsFactors = FALSE)
  q <- q[round(q$prob, 4) %in% LEVELS, ]

  rows <- list()
  for (case in unique(q$case)) for (N in sort(unique(q$N))) for (test in unique(q$test)) for (lv in LEVELS) {
    sub <- q[q$case == case & q$N == N & q$test == test & abs(q$prob - lv) < 1e-9, ]
    sub <- sub[order(sub$T), ]
    if (nrow(sub) < 8) next
    s <- fit_surface(sub$T, sub$cv, sub$mcse)
    if (!(is.finite(s$beta[1]) && is.finite(s$se[1]))) next
    tail <- sub[sub$T >= TAIL_T, ]
    if (nrow(tail) < 1) next
    resid <- (tail$cv - vapply(tail$T, function(TT) cv_at(s, TT), numeric(1))) / tail$mcse
    last <- sub[nrow(sub), ]
    gap <- s$beta[1] - last$cv
    se_gap <- sqrt(s$se[1]^2 + last$mcse^2)
    A <- max(abs(resid)) < 2
    B <- abs(gap) <= 2 * se_gap
    rows[[length(rows) + 1]] <- data.frame(
      case = case, N = N, test = test, level = lv,
      beta_inf = s$beta[1], se_inf = s$se[1], T_max = last$T,
      cv_at_Tmax = last$cv, mcse_Tmax = last$mcse, gap = gap,
      gap_in_se = abs(gap) / se_gap, max_tail_resid = max(abs(resid)),
      n_tail = nrow(tail), sigma = s$sigma, passes_A = A, passes_B = B
    )
  }
  if (length(rows) == 0) stop("no surfaces with at least 8 sample sizes -- is the grid finished?")
  out <- do.call(rbind, rows)
  write.csv(out, file.path(TABLES, "tail_fit_check.csv"), row.names = FALSE)

  fa <- mean(out$passes_A); fb <- mean(out$passes_B)
  cat(sprintf("\nDoes the grid reach the limit?  (tail from T >= %d)\n", TAIL_T))
  cat(strrep("=", 64), "\n")
  cat(sprintf("  combinations                          : %d\n", nrow(out)))
  cat(sprintf("  A: max |tail standardized residual|   : median %.2f, share passing %.2f\n",
              stats::median(out$max_tail_resid), fa))
  cat(sprintf("  B: |beta_inf - c(T_max)| in SE        : median %.2f, share passing %.2f\n",
              stats::median(out$gap_in_se), fb))
  cat(sprintf("  median |gap| in critical-value units  : %.3f\n", stats::median(abs(out$gap))))
  cat(strrep("=", 64), "\n")
  verdict <- if (fa > 0.5 && fb > 0.5) "the grid reaches the limit; beta_inf is interpolated, not extrapolated."
    else if (fa > 0.5) "surface fits the tail but beta_inf sits beyond the data -- extend the grid."
    else if (fb > 0.5) "beta_inf is close to the data but the tail is not on the surface -- respecify."
    else "neither criterion met; the column has not settled on this grid."
  cat(sprintf("  VERDICT: %s\n", verdict))
  cat("\n  worst extrapolation gaps:\n")
  ord <- order(-out$gap_in_se)
  for (i in head(ord, 5)) {
    r <- out[i, ]
    cat(sprintf("    case %-2s N=%d %-4s %2.0f%%  beta_inf %8.3f vs c(%d) %8.3f  gap %+.3f (%.2f se)\n",
                r$case, r$N, r$test, 100 * r$level, r$beta_inf, r$T_max, r$cv_at_Tmax, r$gap, r$gap_in_se))
  }
  cat("\n  wrote tail_fit_check.csv\n")
}

main()
