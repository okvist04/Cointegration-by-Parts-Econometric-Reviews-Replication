# =====================================================================
# verify_engine.R  --  run this before any production run.
#
#   Rscript verify_engine.R
#
# Checks, in order:
#   1. fast engine vs naive reference, window by window, all three
#      deterministic cases, N = 1 and N = 3.
#   3. the three scan statistics vs brute force over the same window set.
#   4. how much the common-sample BIC fix moves the statistic (informational).
#   5. timing of the fast engine, printed as s/rep.
#
# (Numbering matches verify_engine.jl; step 2 there compares against a
# legacy Julia file this repository doesn't have, so it's omitted here.
# The correctness gate this repository actually needs -- fast engine
# agrees with a naive brute-force reference to floating-point precision
# -- is checks 1 and 3.)
# =====================================================================

source("R/engine.R")

fails <- 0
check <- function(name, ok) {
  if (!isTRUE(ok)) fails <<- fails + 1
  cat(sprintf("  %-58s %s\n", name, if (isTRUE(ok)) "ok" else "FAIL"))
}

cat("\n1. fast engine vs naive reference\n")
for (case in c("n", "c", "ct")) {
  for (N in c(1, 3)) {
    set.seed(11)
    T <- 120
    d <- dgp_null(T, N)
    worst <- 0; nbad <- 0; ncmp <- 0
    for (s in seq(1, 60, by = 8)) {
      for (e in seq(s + 24, T, by = 9)) {
        rf <- window_stat_cpp(d$y, d$x, case_to_code(case), s, e)
        rn <- window_stat_naive_cpp(d$y, d$x, case_to_code(case), s, e, common_sample = TRUE)
        ncmp <- ncmp + 1
        if (!(is.finite(rf$tstat) && is.finite(rn$tstat))) { nbad <- nbad + 1; next }
        relerr <- abs(rf$tstat - rn$tstat) / max(1.0, abs(rn$tstat))
        worst <- max(worst, relerr)
        if (rf$lag != rn$lag) nbad <- nbad + 1
      }
    }
    check(sprintf("case %s, N=%d (%d windows, max rel err %.2g)", case, N, ncmp, worst),
          worst < 1e-8 && nbad == 0)
  }
}

cat("\n3. scan statistics vs brute force\n")
{
  set.seed(7)
  T <- 90; N <- 1
  d <- dgp_null(T, N)
  m <- min_window(T, 0.15, 20)
  r <- scan_all3_data_cpp(d$y, d$x, case_to_code("c"), m)
  g <- Inf; f <- Inf; b <- Inf
  for (s in 1:(T - m + 1)) {
    for (e in (s + m - 1):T) {
      w <- window_stat_cpp(d$y, d$x, case_to_code("c"), s, e)
      if (!is.finite(w$tstat)) next
      if (w$tstat < g) g <- w$tstat
      if (s == 1 && w$tstat < f) f <- w$tstat
      if (e == T && w$tstat < b) b <- w$tstat
    }
  }
  check("GIEG", abs(r$gieg - g) < 1e-9)
  check("FIEG", abs(r$fieg - f) < 1e-9)
  check("BIEG", abs(r$bieg - b) < 1e-9)
  check("nesting GIEG <= min(FIEG, BIEG)", r$gieg <= min(r$fieg, r$bieg) + 1e-9)
}

cat("\n4. effect of the common-sample BIC fix (informational, not a pass/fail check)\n")
{
  set.seed(21)
  T <- 200
  N <- 1
  diffs <- c(); dp <- 0
  for (rep in 1:20) {
    d <- dgp_null(T, N)
    for (win in list(c(1, 60), c(1, 200), c(40, 120), c(80, 200))) {
      a <- window_stat_naive_cpp(d$y, d$x, case_to_code("c"), win[1], win[2], common_sample = TRUE)
      b <- window_stat_naive_cpp(d$y, d$x, case_to_code("c"), win[1], win[2], common_sample = FALSE)
      if (!(is.finite(a$tstat) && is.finite(b$tstat))) next
      diffs <- c(diffs, a$tstat - b$tstat)
      if (a$lag != b$lag) dp <- dp + 1
    }
  }
  cat(sprintf("  mean diff %+.4f, max |diff| %.4f, lag order differs in %d of %d windows\n",
              mean(diffs), max(abs(diffs)), dp, length(diffs)))
}

cat("\n5. timing of the fast engine (single core, s/rep)\n")
{
  for (T in c(100, 200, 400, 800)) {
    N <- 1
    m <- min_window(T, 0.15, 20)
    set.seed(1)
    d <- dgp_null(T, N)  # warm-up
    invisible(scan_all3_data_cpp(d$y, d$x, case_to_code("c"), m))
    nrep <- if (T <= 200) 20 else if (T <= 400) 8 else 3
    t0 <- Sys.time()
    for (i in seq_len(nrep)) {
      d <- dgp_null(T, N)
      invisible(scan_all3_data_cpp(d$y, d$x, case_to_code("c"), m))
    }
    el <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("  T=%4d  %8.4f s/rep\n", T, el / nrep))
  }
}

cat(if (fails == 0) "\nALL CHECKS PASSED\n" else sprintf("\n%d CHECK(S) FAILED\n", fails))
quit(status = if (fails == 0) 0 else 1)
