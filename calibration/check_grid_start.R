# =====================================================================
# check_grid_start.R  --  is the smallest T pulling the response surface?
#
#   Rscript check_grid_start.R
#
# The smallest sample size in the grid is simultaneously the noisiest
# row and the row doing most of the work identifying finite-sample
# curvature (recovered from curvature in 1/T). This checks whether it's
# consistent with the rest of the surface.
#
# DECISION RULE (fixed in advance):
#   A. Refit every (case, N, test, level) surface with DROP_T dropped.
#      If |beta_inf(full) - beta_inf(dropped)| exceeds one SE of
#      beta_inf(full) in MORE THAN HALF the combinations, move the grid
#      start up and say so in the caption.
#   B. Independently, if DROP_T's standardized residual in the full fit
#      exceeds 2.5 in more than half the combinations, drop it regardless.
#   Otherwise keep it; report both fits as a robustness note.
#
# Environment: TABLES [results/tables]   DROP_T [100]
# =====================================================================

source("R/response_surface.R")

TABLES <- Sys.getenv("TABLES", "results/tables")
DROP_T <- as.integer(Sys.getenv("DROP_T", "100"))
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
    if (!(DROP_T %in% sub$T)) next
    full <- fit_surface(sub$T, sub$cv, sub$mcse)
    keep_rows <- sub[sub$T != DROP_T, ]
    dropfit <- fit_surface(keep_rows$T, keep_rows$cv, keep_rows$mcse)
    if (!(is.finite(full$beta[1]) && is.finite(dropfit$beta[1]))) next
    row <- sub[sub$T == DROP_T, ][1, ]
    r_std <- (row$cv - cv_at(full, row$T)) / row$mcse
    shift <- full$beta[1] - dropfit$beta[1]
    rows[[length(rows) + 1]] <- data.frame(
      case = case, N = N, test = test, level = lv,
      beta_inf_full = full$beta[1], se_full = full$se[1],
      beta_inf_drop = dropfit$beta[1], shift = shift,
      shift_in_se = abs(shift) / full$se[1], resid_T_std = r_std
    )
  }
  if (length(rows) == 0) stop("no comparable surfaces -- is the grid finished?")
  out <- do.call(rbind, rows)
  write.csv(out, file.path(TABLES, "grid_start_check.csv"), row.names = FALSE)

  frac_a <- mean(out$shift_in_se > 1)
  frac_b <- mean(abs(out$resid_T_std) > 2.5)
  cat(sprintf("\nT = %d sensitivity of the response surface\n", DROP_T))
  cat(strrep("=", 62), "\n")
  cat(sprintf("  combinations compared                : %d\n", nrow(out)))
  cat(sprintf("  median |shift| in standard errors    : %.2f\n", stats::median(out$shift_in_se)))
  cat(sprintf("  share with |shift| > 1 se            : %.2f   (rule A: >0.5 drops it)\n", frac_a))
  cat(sprintf("  median |standardized residual| at T  : %.2f\n", stats::median(abs(out$resid_T_std))))
  cat(sprintf("  share with |residual| > 2.5          : %.2f   (rule B: >0.5 drops it)\n", frac_b))
  cat(strrep("=", 62), "\n")
  verdict <- (frac_a > 0.5) || (frac_b > 0.5)
  cat(if (verdict) sprintf("  VERDICT: move the grid start above T = %d and say so in the caption.\n", DROP_T)
      else sprintf("  VERDICT: keep T = %d. Report both fits as a robustness note.\n", DROP_T))
  cat("\n  largest shifts:\n")
  ord <- order(-out$shift_in_se)
  for (i in head(ord, 5)) {
    r <- out[i, ]
    cat(sprintf("    case %-2s N=%d %-4s %2.0f%%  beta_inf %8.3f -> %8.3f  (%.2f se)\n",
                r$case, r$N, r$test, 100 * r$level, r$beta_inf_full, r$beta_inf_drop, r$shift_in_se))
  }
  cat("\n  wrote grid_start_check.csv\n")
}

main()
