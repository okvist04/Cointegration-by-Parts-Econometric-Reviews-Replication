# ============================================================
# cointegration_by_parts.R
# Full R translation of the Julia implementation
#
# CHANGES FROM ORIGINAL:
#  - BUG FIX: n_regressors (and beta) were silently ignored in the
#    Monte Carlo null path -- every past run used n_regressors = 1
#    no matter what was requested. Fixed and threaded through properly.
#  - The expensive per-window scan (scan_expanding_windows /
#    scan_all_subsamples, called via residual_df_window /
#    adf_noconst_fast) now has a fast C++ backend in
#    adf_scan_fast.cpp, used by default for the null-calibration path.
#    The original pure-R implementations are KEPT BELOW UNCHANGED so
#    you can validate the fast path against them (see
#    validate_fast_vs_r.R) -- they are simply no longer on the hot path.
#  - calibrate_null_critical_values_parallel() now loops over a
#    n_regressors grid too, adds a bootstrap SE per critical value, and
#    checkpoints (and can resume from) each (T, n_regressors, scan_type)
#    combination independently.
# ============================================================

library(parallel)
library(data.table)
library(pbmcapply)
library(Rcpp)

# Compile the fast scan engine once per session. Path is relative to
# the working directory -- adjust if you keep the .cpp file elsewhere.
sourceCpp("Claude/adf_scan_fast.cpp")

# ============================================================
# simulate_pair_break
# ============================================================

simulate_pair_break <- function(
    T,
    n_regressors = 1,
    frac_coint,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_e = 1,
    rho_stationary = 0.3,
    regime = "coint_to_noncoint",
    seed = NULL
) {

  if (T < 10)
    stop("T must be at least 10.")

  if (frac_coint < 0 || frac_coint > 1)
    stop("frac_coint must be in [0,1].")

  if (abs(rho_stationary) >= 1)
    stop("abs(rho_stationary) must be < 1.")

  if (!(regime %in% c("coint_to_noncoint", "noncoint_to_coint")))
    stop("Invalid regime.")

  if (!is.null(seed))
    set.seed(seed)

  tau <- floor(frac_coint * T)

  X <- matrix(NA_real_, T, n_regressors)

  for (j in seq_len(n_regressors)) {
    X[, j] <- cumsum(sigma_x * rnorm(T))
  }

  if (is.null(beta))
    beta <- rep(1, n_regressors)

  if (length(beta) != n_regressors)
    stop("beta must have length n_regressors.")

  eta <- sigma_e * rnorm(T)

  e <- numeric(T)

  for (t in 2:T) {

    if (regime == "coint_to_noncoint") {

      if (t <= tau) {
        e[t] <- rho_stationary * e[t - 1] + eta[t]
      } else {
        e[t] <- e[t - 1] + eta[t]
      }

    } else {

      if (t <= tau) {
        e[t] <- e[t - 1] + eta[t]
      } else {
        e[t] <- rho_stationary * e[t - 1] + eta[t]
      }
    }
  }

  X <- as.matrix(X)

  if (ncol(X) != length(beta))
    stop("length(beta) must equal ncol(X).")

  y <- as.numeric(alpha + X %*% beta + e)

  list(
    y = y,
    X = X,
    e = e,
    tau = tau,
    frac_coint = frac_coint,
    regime = regime
  )
}

# ============================================================
# simulate_window_break
#
# BUG FIX: original used a free variable `n_regressors` that was not a
# parameter of this function (would error if ever called). Added as a
# proper parameter, default 1, same as simulate_pair_break.
# ============================================================

simulate_window_break <- function(
    T,
    break_fracs,
    n_regressors = 1,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_e = 1,
    rho_stationary = 0.3,
    starting_regime = "coint",
    seed = NULL
) {

  if (T < 10)
    stop("T must be at least 10.")

  if (length(break_fracs) == 0)
    stop("break_fracs must be non-empty.")

  if (any(break_fracs <= 0 | break_fracs >= 1))
    stop("All break_fracs must be in (0,1).")

  if (is.unsorted(break_fracs))
    stop("break_fracs must be sorted.")

  if (!(starting_regime %in% c("coint", "noncoint")))
    stop("Invalid starting_regime.")

  if (!is.null(seed))
    set.seed(seed)

  break_idx <- floor(break_fracs * T)

  X <- matrix(NA_real_, T, n_regressors)

  for (j in seq_len(n_regressors)) {
    X[, j] <- cumsum(sigma_x * rnorm(T))
  }

  if (is.null(beta))
    beta <- rep(1, n_regressors)

  if (length(beta) != n_regressors)
    stop("beta must have length n_regressors.")

  eta <- sigma_e * rnorm(T)

  e <- numeric(T)

  for (t in 2:T) {

    seg <- sum(break_idx < t)

    is_coint <- if (starting_regime == "coint") {
      seg %% 2 == 0
    } else {
      seg %% 2 == 1
    }

    if (is_coint) {
      e[t] <- rho_stationary * e[t - 1] + eta[t]
    } else {
      e[t] <- e[t - 1] + eta[t]
    }
  }

  X <- as.matrix(X)

  if (ncol(X) != length(beta))
    stop("length(beta) must equal ncol(X).")

  y <- as.numeric(alpha + X %*% beta + e)

  list(
    y = y,
    X = X,
    e = e,
    break_idx = break_idx,
    break_fracs = break_fracs,
    starting_regime = starting_regime
  )
}

# ============================================================
# Fast ADF without deterministic terms (PURE R -- kept for validation)
# ============================================================

adf_noconst_fast <- function(
    series,
    max_lag = NULL,
    criterion = "bic"
) {

  T <- length(series)

  if (T < 10) {
    return(list(tstat = NaN, lag = -1))
  }

  if (is.null(max_lag)) {
    max_lag <- floor(12 * (T / 100)^(1 / 4))
  }

  max_lag <- max(0, max_lag)

  dy <- diff(series)
  ylag <- series[-length(series)]

  best_ic <- Inf
  best_t <- NaN
  best_p <- -1

  for (p in 0:max_lag) {

    n <- length(dy) - p

    if (n <= 0)
      next

    Y <- dy[(p + 1):length(dy)]
    yreg <- ylag[(p + 1):length(ylag)]

    if (p == 0) {

      X <- matrix(yreg, ncol = 1)

    } else {

      Xlags <- matrix(NA_real_, nrow = n, ncol = p)

      for (j in 1:p) {
        Xlags[, j] <- dy[(p + 1 - j):(length(dy) - j)]
      }

      X <- cbind(yreg, Xlags)
    }

    k <- ncol(X)

    if (qr(X)$rank < k)
      next

    if (n <= k)
      next

    fit <- lm.fit(X, Y)

    beta_hat <- fit$coefficients
    resid <- fit$residuals

    rss <- sum(resid^2)

    if (rss <= 0)
      next

    sigma2 <- rss / (n - k)

    xtx_inv <- solve(crossprod(X))

    se_phi <- sqrt(sigma2 * xtx_inv[1, 1])

    if (se_phi <= 0)
      next

    tstat <- beta_hat[1] / se_phi

    ic <- if (criterion == "aic") {
      n * log(rss / n) + 2 * k
    } else if (criterion == "bic") {
      n * log(rss / n) + log(n) * k
    } else {
      stop("criterion must be 'aic' or 'bic'")
    }

    if (ic < best_ic) {
      best_ic <- ic
      best_t <- tstat
      best_p <- p
    }
  }

  list(
    tstat = best_t,
    lag = best_p
  )
}

# ============================================================
# residual_df_window (PURE R -- kept for validation)
# ============================================================

residual_df_window <- function(
    y_sub,
    X_sub,
    include_intercept = TRUE,
    max_lag = NULL
) {

  ## Convert to matrix
  X <- as.matrix(X_sub)

  if (nrow(X) != length(y_sub))
    stop("y_sub and X_sub must have the same number of observations.")

  if (length(y_sub) < 10)
    return(list(tstat = NaN, lag = -1))

  if (include_intercept) {
    X <- cbind(1, X)
  }

  if (qr(X)$rank < ncol(X))
    return(list(tstat = NaN, lag = -1))

  fit <- lm.fit(X, y_sub)

  e_hat <- fit$residuals

  adf_noconst_fast(
    e_hat,
    max_lag = max_lag,
    criterion = "bic"
  )
}

# ============================================================
# argmin_finite
#
# Index of the minimum value in v, ignoring non-finite entries (NA,
# NaN, Inf). Returns 0 if no finite entries exist, so callers can
# check `idx == 0` rather than handle an error -- convenient inside
# the scan loops, where a window can legitimately produce no finite
# statistic (e.g. rank-deficient regression) and the loop should just
# skip it rather than stop.
#
# Args:
#   v -- numeric vector, possibly containing NA/NaN/Inf.
# Returns:
#   Integer index of the minimum finite value, or 0 if none exist.
# ============================================================

argmin_finite <- function(v) {

  idx <- which(is.finite(v))

  if (length(idx) == 0)
    return(0)

  idx[which.min(v[idx])]
}

# ============================================================
# scan_expanding_windows (PURE R -- kept for validation)
# ============================================================

scan_expanding_windows <- function(
    y,
    X,
    direction = "forward",
    min_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL
) {

  X <- as.matrix(X)

  if (nrow(X) != length(y))
    stop("y and x must have the same number of observations.")

  if (!(direction %in% c("forward", "backward")))
    stop("direction must be forward/backward")

  T <- length(y)

  if (min_obs < 10 || min_obs > T)
    stop("Invalid min_obs")

  n_windows <- T - min_obs + 1

  start_idx <- integer(n_windows)
  end_idx <- integer(n_windows)
  n_obs <- integer(n_windows)
  tstats <- numeric(n_windows)
  lags <- integer(n_windows)

  for (i in 1:n_windows) {

    if (direction == "forward") {

      s <- 1
      e <- min_obs + i - 1

    } else {

      s <- T - (min_obs + i - 1) + 1
      e <- T
    }

    adf_res <- residual_df_window(
      y[s:e],
      X[s:e, , drop = FALSE],
      include_intercept = include_intercept,
      max_lag = max_lag
    )

    start_idx[i] <- s
    end_idx[i] <- e
    n_obs[i] <- e - s + 1
    tstats[i] <- adf_res$tstat
    lags[i] <- adf_res$lag
  }

  min_id <- argmin_finite(tstats)

  list(
    scan = data.table(
      window_id = 1:n_windows,
      start = start_idx,
      stop = end_idx,
      n_obs = n_obs,
      tstat = tstats,
      lag = lags
    ),
    min_t = ifelse(min_id == 0, NaN, tstats[min_id]),
    min_window_id = min_id,
    min_start = ifelse(min_id == 0, 0, start_idx[min_id]),
    min_stop = ifelse(min_id == 0, 0, end_idx[min_id])
  )
}

# ============================================================
# scan_all_subsamples (PURE R -- kept for validation)
# ============================================================

scan_all_subsamples <- function(
    y,
    X,
    min_start = 1,
    min_window_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL
) {

  X <- as.matrix(X)

  if (nrow(X) != length(y))
    stop("y and x must have the same number of observations.")

  T <- length(y)

  s_max <- T - min_window_obs + 1

  total_windows <- 0

  for (s in min_start:s_max) {
    total_windows <- total_windows +
      (T - (s + min_window_obs - 1) + 1)
  }

  window_id <- integer(total_windows)
  start_idx <- integer(total_windows)
  end_idx <- integer(total_windows)
  n_obs <- integer(total_windows)
  tstats <- numeric(total_windows)
  lags <- integer(total_windows)

  k <- 0

  for (s in min_start:s_max) {

    e_min <- s + min_window_obs - 1

    for (e in e_min:T) {

      k <- k + 1

      adf_res <- residual_df_window(
        y[s:e],
        X[s:e, , drop = FALSE],
        include_intercept = include_intercept,
        max_lag = max_lag
      )

      window_id[k] <- k
      start_idx[k] <- s
      end_idx[k] <- e
      n_obs[k] <- e - s + 1
      tstats[k] <- adf_res$tstat
      lags[k] <- adf_res$lag
    }
  }

  min_id <- argmin_finite(tstats)

  list(
    scan = data.table(
      window_id = window_id,
      start = start_idx,
      stop = end_idx,
      n_obs = n_obs,
      tstat = tstats,
      lag = lags
    ),
    min_t = ifelse(min_id == 0, NaN, tstats[min_id]),
    min_window_id = min_id,
    min_start = ifelse(min_id == 0, 0, start_idx[min_id]),
    min_stop = ifelse(min_id == 0, 0, end_idx[min_id])
  )
}

# ============================================================
# simulate_pair_null
# ============================================================

# ============================================================
# simulate_pair_null
#
# Simulates one draw of a NON-cointegrated pair (y, X) under the null
# hypothesis: X is a matrix of independent random walks, and y is
# alpha + X %*% beta + u, where u is an INDEPENDENT random walk (not a
# stationary combination of X). Because u itself is I(1), y and X are
# NOT cointegrated by construction -- this is exactly the null
# scenario the calibration grid (calibrate_null_critical_values_parallel)
# needs to simulate from repeatedly to build the distribution of the
# minimum test statistic under "no cointegration".
#
# Args:
#   T              -- sample size (>= 10).
#   n_regressors   -- number of columns in X.
#   alpha          -- intercept in the y equation.
#   beta           -- coefficient vector on X (length n_regressors);
#                     defaults to a vector of 1s if NULL.
#   sigma_x        -- innovation SD for each random walk in X.
#   sigma_u        -- innovation SD for the independent random walk u.
#   seed           -- optional seed for reproducibility.
# Returns:
#   A list with y (numeric vector), X (T x n_regressors matrix), and
#   u (the independent random-walk error, for diagnostics).
# ============================================================

simulate_pair_null <- function(
    T,
    n_regressors = 1,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_u = 1,
    seed = NULL
) {

  if (T < 10)
    stop("T must be at least 10.")

  if (!is.null(seed))
    set.seed(seed)

  X <- matrix(NA_real_, T, n_regressors)

  for (j in seq_len(n_regressors)) {
    X[, j] <- cumsum(sigma_x * rnorm(T))
  }

  if (is.null(beta))
    beta <- rep(1, n_regressors)

  if (length(beta) != n_regressors)
    stop("beta must have length n_regressors.")

  u <- cumsum(sigma_u * rnorm(T))

  X <- as.matrix(X)

  if (ncol(X) != length(beta))
    stop("length(beta) must equal ncol(X).")

  y <- as.numeric(alpha + X %*% beta + u)

  list(
    y = y,
    X = X,
    u = u
  )
}

# ============================================================
# one null replication -- FAST (C++-backed) path, used by default
#
# BUG FIX: n_regressors and beta are now actually threaded through to
# simulate_pair_null() (previously hardcoded to n_regressors = 1,
# beta = NULL regardless of what was requested).
# ============================================================

minstat_one_rep_null <- function(
    T,
    n_regressors = 1,
    scan_type = "fw",
    min_obs = 40,
    min_start = 1,
    min_window_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_u = 1,
    seed = NULL
) {

  d <- simulate_pair_null(
    T = T,
    n_regressors = n_regressors,
    alpha = alpha,
    beta = beta,
    sigma_x = sigma_x,
    sigma_u = sigma_u,
    seed = seed
  )

  ml <- if (is.null(max_lag)) -1L else as.integer(max_lag)

  minstat_scan_cpp(
    y_r = d$y,
    X_r = as.matrix(d$X),
    scan_type = scan_type,
    min_obs = as.integer(min_obs),
    min_start = as.integer(min_start),
    min_window_obs = as.integer(min_window_obs),
    include_intercept = include_intercept,
    max_lag = ml
  )
}

# ============================================================
# one null replication -- PURE R reference path (kept for validation
# against the fast path; same bug fix applied so the comparison is
# apples-to-apples)
# ============================================================

minstat_one_rep_null_r <- function(
    T,
    n_regressors = 1,
    scan_type = "fw",
    min_obs = 40,
    min_start = 1,
    min_window_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_u = 1,
    seed = NULL
) {

  d <- simulate_pair_null(
    T = T,
    n_regressors = n_regressors,
    alpha = alpha,
    beta = beta,
    sigma_x = sigma_x,
    sigma_u = sigma_u,
    seed = seed
  )

  if (scan_type == "fw") {

    res <- scan_expanding_windows(
      d$y, d$X,
      direction = "forward",
      min_obs = min_obs,
      include_intercept = include_intercept,
      max_lag = max_lag
    )

  } else if (scan_type == "bw") {

    res <- scan_expanding_windows(
      d$y, d$X,
      direction = "backward",
      min_obs = min_obs,
      include_intercept = include_intercept,
      max_lag = max_lag
    )

  } else {

    res <- scan_all_subsamples(
      d$y, d$X,
      min_start = min_start,
      min_window_obs = min_window_obs,
      include_intercept = include_intercept,
      max_lag = max_lag
    )
  }

  res$min_t
}

# ============================================================
# simulate_minstats_null
#
# BUG FIX: n_regressors and beta now actually reach the simulator.
# Added `engine` so you can switch to the pure-R path for validation.
# ============================================================

simulate_minstats_null <- function(
    T,
    n_rep,
    n_regressors = 1,
    scan_type = "fw",
    min_obs = 40,
    min_start = 1,
    min_window_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_u = 1,
    base_seed = 20260429,
    n_cores = max(1, detectCores() - 1),
    engine = c("cpp", "r")
) {

  engine <- match.arg(engine)
  worker <- if (engine == "cpp") minstat_one_rep_null else minstat_one_rep_null_r

  seeds <- base_seed + seq_len(n_rep)

  res <- pbmclapply(
    seeds,
    mc.cores = n_cores,
    mc.preschedule = FALSE, # task cost varies a lot with T / scan_type;
                             # dynamic scheduling balances cores better
    FUN = function(s) {

      worker(
        T = T,
        n_regressors = n_regressors,
        scan_type = scan_type,
        min_obs = min_obs,
        min_start = min_start,
        min_window_obs = min_window_obs,
        include_intercept = include_intercept,
        max_lag = max_lag,
        alpha = alpha,
        beta = beta,
        sigma_x = sigma_x,
        sigma_u = sigma_u,
        seed = s
      )
    }
  )

  unlist(res)
}

# ============================================================
# empirical critical values
# ============================================================

# ============================================================
# empirical_critical_values
#
# Given a vector of simulated minimum test statistics under the null
# (one per replication, from simulate_minstats_null()), returns the
# empirical quantiles at the requested probabilities -- these are the
# raw (pre-response-surface) critical values for one specific
# (T, n_regressors, scan_type) combination. fit_response_surfaces()
# later smooths these across the T grid via a response-surface
# regression, but this function is what produces the per-T inputs to
# that regression in the first place.
#
# Args:
#   minstats -- numeric vector of simulated minimum statistics (one
#               per replication); non-finite entries are dropped.
#   probs    -- significance levels to compute quantiles at.
# Returns:
#   A list with probs, crit_vals (the corresponding quantiles),
#   n_finite (usable replications), and n_total (all replications,
#   for diagnosing how many were dropped as non-finite).
# ============================================================

empirical_critical_values <- function(
    minstats,
    probs = c(0.01, 0.05, 0.10)
) {

  finite_stats <- minstats[is.finite(minstats)]

  if (length(finite_stats) == 0)
    stop("No finite statistics available.")

  qvals <- quantile(
    finite_stats,
    probs = probs,
    names = FALSE
  )

  list(
    probs = probs,
    crit_vals = qvals,
    n_finite = length(finite_stats),
    n_total = length(minstats)
  )
}

# ============================================================
# bootstrap SE for the empirical critical values
#
# Nonparametric bootstrap: resample the finite min-stats with
# replacement n_boot times and take the SD of the resulting quantile
# estimates. Cheap relative to the simulation itself (n_boot resamples
# of a length-n_rep vector), so it's done by default.
# ============================================================

bootstrap_quantile_se <- function(
    minstats,
    probs = c(0.01, 0.05, 0.10),
    n_boot = 500,
    seed = NULL
) {

  x <- minstats[is.finite(minstats)]
  n <- length(x)

  if (n == 0)
    return(rep(NA_real_, length(probs)))

  if (!is.null(seed))
    set.seed(seed)

  boot_q <- matrix(NA_real_, nrow = n_boot, ncol = length(probs))

  for (b in seq_len(n_boot)) {
    xb <- x[sample.int(n, n, replace = TRUE)]
    boot_q[b, ] <- quantile(xb, probs = probs, names = FALSE)
  }

  apply(boot_q, 2, sd)
}

# ============================================================
# empirical p-value
# ============================================================

# ============================================================
# empirical_pvalue
#
# One-sided empirical p-value for an observed test statistic against
# the simulated null distribution of minimum statistics: the
# proportion of null replications at least as extreme (<=, since the
# test rejects for very negative statistics) as test_stat. Uses the
# standard (count + 1) / (n + 1) correction so the p-value is never
# exactly 0 even if test_stat is more extreme than every simulated
# replication.
#
# Args:
#   test_stat -- the observed test statistic to evaluate.
#   minstats  -- numeric vector of simulated null minimum statistics.
# Returns:
#   A single numeric p-value in (0, 1].
# ============================================================

empirical_pvalue <- function(
    test_stat,
    minstats
) {

  finite_stats <- minstats[is.finite(minstats)]

  if (length(finite_stats) == 0)
    stop("No finite statistics available.")

  (sum(finite_stats <= test_stat) + 1) /
    (length(finite_stats) + 1)
}

# ============================================================
# calibrate_null_critical_values_parallel
#
# Now loops over T_grid x n_regressors_grid x scan_types. Each
# (T, n_regressors, scan_type) combination is checkpointed to its own
# file; if that file already exists (e.g. from a previous, interrupted
# run) it's loaded instead of recomputed, so the whole grid is safe to
# re-launch after a crash/timeout without losing completed work.
# ============================================================

calibrate_null_critical_values_parallel <- function(
    T_grid = c(200, 500, 1000, 2000, 5000, 10000),
    n_rep = 10000,
    n_regressors_grid = 1,
    scan_types = c("fw", "bw", "all"),
    probs = c(0.01, 0.05, 0.10),
    min_obs = 40,
    min_start = 1,
    min_window_obs = 40,
    include_intercept = TRUE,
    max_lag = NULL,
    alpha = 0,
    beta = NULL,
    sigma_x = 1,
    sigma_u = 1,
    base_seed = 20260429,
    n_cores = max(1, detectCores() - 1),
    n_boot_se = 500,
    checkpoint_dir = ".",
    engine = c("cpp", "r")
) {

  engine <- match.arg(engine)

  if (!dir.exists(checkpoint_dir))
    dir.create(checkpoint_dir, recursive = TRUE)

  combos <- expand.grid(
    T = T_grid,
    n_regressors = n_regressors_grid,
    scan_type = scan_types,
    stringsAsFactors = FALSE
  )

  results <- vector("list", nrow(combos))

  pb <- txtProgressBar(min = 0, max = nrow(combos), style = 3)

  # Seeds are offset deterministically per combo (not just per T) so
  # that adding n_regressors to the grid doesn't shift the seeds used
  # by combinations that were already run and checkpointed.
  seed_for_combo <- function(i) {
    base_seed + (i - 1) * n_rep
  }

  for (i in seq_len(nrow(combos))) {

    T_i <- combos$T[i]
    k_i <- combos$n_regressors[i]
    st_i <- combos$scan_type[i]

    ckpt_file <- file.path(
      checkpoint_dir,
      sprintf("checkpoint_T%d_k%d_%s.csv", T_i, k_i, st_i)
    )

    if (file.exists(ckpt_file)) {

      cat(sprintf(
        "\n[skip, already done] T=%d n_regressors=%d scan_type=%s -> %s\n",
        T_i, k_i, st_i, ckpt_file
      ))

      results[[i]] <- fread(ckpt_file)

    } else {

      cat(sprintf(
        "\n[running] T=%d n_regressors=%d scan_type=%s\n",
        T_i, k_i, st_i
      ))

      minstats <- simulate_minstats_null(
        T = T_i,
        n_rep = n_rep,
        n_regressors = k_i,
        scan_type = st_i,
        min_obs = min_obs,
        min_start = min_start,
        min_window_obs = min_window_obs,
        include_intercept = include_intercept,
        max_lag = max_lag,
        alpha = alpha,
        beta = beta,
        sigma_x = sigma_x,
        sigma_u = sigma_u,
        base_seed = seed_for_combo(i),
        n_cores = n_cores,
        engine = engine
      )

      cv <- empirical_critical_values(minstats, probs = probs)

      se <- bootstrap_quantile_se(
        minstats,
        probs = probs,
        n_boot = n_boot_se,
        seed = seed_for_combo(i) + 1L
      )

      combo_result <- data.table(
        T = T_i,
        n_regressors = k_i,
        scan_type = st_i,
        prob = cv$probs,
        crit_val = cv$crit_vals,
        crit_val_se = se,
        n_finite = cv$n_finite,
        n_total = cv$n_total
      )

      fwrite(combo_result, ckpt_file)

      results[[i]] <- combo_result
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  rbindlist(results)
}
