# =====================================================================
# paper/ndf_check.R  --  diagnostics for the numerical distribution
# function: the probability grid, on-grid accuracy, off-grid
# calibration, and the beta_inf extrapolation check for GIEG.
#
#   Rscript paper/ndf_check.R
# =====================================================================

source("R/response_surface.R")

SIM <- "."
q   <- read.csv(file.path(SIM, "results", "tables", "quantiles.csv"), stringsAsFactors = FALSE)
ndf <- read.csv(file.path(SIM, "results", "tables", "surface_ndf.csv"), stringsAsFactors = FALSE)
TESTS <- c("FIEG", "BIEG", "GIEG")

# ---- A. the probability grid --------------------------------------
pg <- sort(unique(ndf$prob))
cat(sprintf("A. NDF grid: %d points, min %.5f, max %.5f\n", length(pg), min(pg), max(pg)))
cat(sprintf("   points <= 0.01: %d ; in (0.01,0.20]: %d ; > 0.20: %d\n",
            sum(pg <= 0.01), sum(pg > 0.01 & pg <= 0.20), sum(pg > 0.20)))
cat("   first 6:", round(head(pg, 6), 5), "\n")

surfaces_for <- function(case, N, tst) {
  s <- ndf[ndf$case == case & ndf$N == N & ndf$test == tst, ]
  s <- s[order(s$prob), ]
  list(probs = s$prob, surfaces = lapply(seq_len(nrow(s)), function(i) {
    r <- s[i, ]
    list(beta = c(r$beta_inf, r$b1, r$b2, r$b3), se = c(r$se_inf, r$se1, r$se2, r$se3),
         sigma = r$sigma, cubic = r$cubic, dof = r$dof, nobs = r$nobs)
  }))
}

# ---- B. on-grid accuracy: p at the simulated quantile --------------
cat("\nB. NDF evaluated at the simulated quantile (target = nominal level)\n")
rowsB <- list()
for (case in c("c", "ct")) for (N in 1:3) for (tst in TESTS) {
  sf <- surfaces_for(case, N, tst)
  if (length(sf$probs) == 0) next
  for (a in c(0.01, 0.05, 0.10)) {
    sub <- q[q$case == case & q$N == N & q$test == tst & abs(q$prob - a) < 1e-9, ]
    sub <- sub[order(sub$T), ]
    d <- vapply(seq_len(nrow(sub)), function(i) {
      abs(pvalue_from_surface(sf$probs, sf$surfaces, sub$T[i], sub$cv[i]) - a)
    }, numeric(1))
    if (length(d) == 0) next
    rowsB[[length(rowsB) + 1]] <- data.frame(case = case, N = N, test = tst, level = a,
                                              maxdev = max(d), meandev = mean(d))
  }
}
rowsB_df <- do.call(rbind, rowsB)
for (a in c(0.01, 0.05, 0.10)) {
  s <- rowsB_df[abs(rowsB_df$level - a) < 1e-9, ]
  cat(sprintf("   level %.2f : mean |dev| %.4f, max |dev| %.4f\n", a, mean(s$meandev), max(s$maxdev)))
}
for (tst in TESTS) {
  s <- rowsB_df[rowsB_df$test == tst, ]
  cat(sprintf("   %-4s      : mean |dev| %.4f, max |dev| %.4f\n", tst, mean(s$meandev), max(s$maxdev)))
}

# ---- C. off-grid: distribution of NDF p-values on stored null draws -
cat("\nC. NDF p-values on null draws at sample sizes off the fitting grid\n")
for (case in c("c", "ct")) for (T in c(90, 140, 350, 600)) for (N in c(1, 3)) {
  f <- file.path(SIM, "results", "sizepower", sprintf("size_case-%s_T-%d_N-%d.rds", case, T, N))
  if (!file.exists(f)) next
  st <- readRDS(f)$stats
  for (tst in TESTS) {
    sf <- surfaces_for(case, N, tst)
    if (length(sf$probs) == 0) next
    v <- st[[tolower(tst)]]
    pv <- vapply(v, function(x) pvalue_from_surface(sf$probs, sf$surfaces, T, x), numeric(1))
    cat(sprintf("   case %-2s T=%4d N=%d %-4s  R=%d  P(p<=.01)=%.3f  P(p<=.05)=%.3f  P(p<=.10)=%.3f\n",
                case, T, N, tst, length(pv), mean(pv <= 0.01, na.rm = TRUE),
                mean(pv <= 0.05, na.rm = TRUE), mean(pv <= 0.10, na.rm = TRUE)))
  }
}

# ---- D. sigma-hat and the beta_inf extrapolation --------------------
cat("\nD. response surface fit at 1/5/10%\n")
main <- ndf[round(ndf$prob, 5) %in% c(0.01, 0.05, 0.10), ]
for (case in c("c", "ct")) for (tst in TESTS) {
  s <- main[main$case == case & main$test == tst, ]
  if (nrow(s) == 0) next
  cat(sprintf("   case %-2s %-4s  sigma median %.2f  range %.2f-%.2f  cubic %d/%d  se(beta_inf) median %.3f max %.3f\n",
              case, tst, stats::median(s$sigma), min(s$sigma), max(s$sigma), sum(s$cubic), nrow(s),
              stats::median(s$se_inf), max(s$se_inf)))
}
cat("\n   beta_inf versus the largest simulated T:\n")
for (case in c("c", "ct")) for (tst in TESTS) for (N in 1:3) {
  s <- main[main$case == case & main$test == tst & main$N == N, ]
  if (nrow(s) == 0) next
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    sub <- q[q$case == case & q$N == N & q$test == tst & abs(q$prob - r$prob) < 1e-9, ]
    Tmax <- max(sub$T); rr <- sub[sub$T == Tmax, ][1, ]
    gap <- r$beta_inf - rr$cv
    cat(sprintf("   case %-2s %-4s N=%d %3d%%  beta_inf %+8.3f  cv(T=%d) %+8.3f  gap %+6.3f  (mcse %.3f, se_inf %.3f)\n",
                case, tst, N, round(100 * r$prob), r$beta_inf, Tmax, rr$cv, gap, rr$mcse, r$se_inf))
  }
}
