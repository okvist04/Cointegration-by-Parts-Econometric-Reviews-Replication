# =====================================================================
# response_surface.R
#
# Port of response_surface.jl: MacKinnon-style response surfaces and
# the numerical distribution function built on top of them.
#
#   c_hat(T) = beta_inf + beta_1/T + beta_2/T^2 + beta_3/T^3
#
# fitted by weighted least squares with weights w_T = se_T^{-2}, se_T
# the bootstrap standard error of the simulated quantile.
#
# NOTE on fidelity: Julia's normcdf()/norminv() are hand-written
# rational approximations (Zelen & Severo / Acklam) used only because
# base Julia has no normal cdf/quantile without an extra package. R's
# built-in pnorm()/qnorm() are exact (to machine precision, via TOMS
# 708 / AS 241) and are used directly here instead -- a strictly more
# accurate substitution, not a simplification. Any difference from the
# Julia numbers is at the level of the two approximations' own error
# (~1e-9), far below anything that matters for a p-value.
# =====================================================================

# Bootstrap standard error of the p-quantile of `draws`. Unlike Julia's
# version (which threads a single rng object through many calls), this
# relies on R's own persistent global RNG stream advancing naturally
# across repeated calls -- call set.seed() ONCE in the calling script
# before looping over cells, exactly as Julia's main() creates its rng
# object once and passes it through every call in the loop.
quantile_se <- function(draws, p, B = 1000) {
  v <- draws[is.finite(draws)]
  n <- length(v)
  if (n <= 10) return(NA_real_)
  qs <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    stats::quantile(v[idx], probs = p, names = FALSE, type = 7)
  }, numeric(1))
  stats::sd(qs)
}

# ---------------------------------------------------------------------
# Weighted least squares fit, mirroring Julia's `wls()`
# ---------------------------------------------------------------------

.wls <- function(X, y, w) {
  sw <- sqrt(w)
  Xw <- X * sw
  yw <- y * sw
  XtX <- crossprod(Xw)
  ch <- tryCatch(chol(XtX), error = function(e) NULL)
  p <- ncol(X)
  if (is.null(ch)) return(list(beta = rep(NA_real_, p), se = rep(NA_real_, p), sigma = NA_real_))
  b <- backsolve(ch, backsolve(ch, crossprod(Xw, yw), transpose = TRUE))
  r <- yw - Xw %*% b
  dof <- length(y) - p
  s2 <- if (dof > 0) sum(r^2) / dof else NA_real_
  V <- chol2inv(ch)
  se <- sqrt(pmax(0, s2 * diag(V)))
  list(beta = as.numeric(b), se = as.numeric(se), sigma = sqrt(s2))
}

# ---------------------------------------------------------------------
# fit_surface: cv_T = beta_inf + beta_1/T + beta_2/T^2 + beta_3/T^3
# ---------------------------------------------------------------------

fit_surface <- function(Ts, cv, se, cubic_t = 1.96) {
  ok <- is.finite(cv) & is.finite(se) & (se > 0)
  T <- as.numeric(Ts[ok]); y <- cv[ok]; w <- se[ok]^(-2)
  n <- length(T)
  empty <- list(beta = rep(NA_real_, 4), se = rep(NA_real_, 4), sigma = NA_real_,
                cubic = FALSE, dof = 0L, nobs = n)
  if (n < 4) return(empty)

  X3 <- cbind(1, 1 / T, 1 / T^2, 1 / T^3)
  f3 <- .wls(X3, y, w)
  if (is.finite(f3$se[4]) && abs(f3$beta[4] / f3$se[4]) >= cubic_t) {
    return(list(beta = f3$beta, se = f3$se, sigma = f3$sigma, cubic = TRUE, dof = n - 4, nobs = n))
  }
  X2 <- X3[, 1:3, drop = FALSE]
  f2 <- .wls(X2, y, w)
  list(beta = c(f2$beta, 0.0), se = c(f2$se, NA_real_), sigma = f2$sigma,
       cubic = FALSE, dof = n - 3, nobs = n)
}

cv_at <- function(surf, T) {
  surf$beta[1] + surf$beta[2] / T + surf$beta[3] / T^2 + surf$beta[4] / T^3
}

# ---------------------------------------------------------------------
# pvalue_from_surface: approximate p-value for `stat` at sample size T,
# from surfaces fitted at the dense probability grid `probs`.
# Interpolation is linear in the normal quantile of p (MacKinnon 1994,
# 1996 style); values off the ends of the grid return the bound.
# ---------------------------------------------------------------------

pvalue_from_surface <- function(probs, surfaces, T, stat) {
  cvals <- vapply(surfaces, cv_at, numeric(1), T = T)
  ok <- is.finite(cvals)
  if (sum(ok) < 4) return(NA_real_)
  c_ok <- cvals[ok]; p_ok <- probs[ok]
  ord <- order(c_ok)
  c_ord <- c_ok[ord]; p_ord <- p_ok[ord]
  if (stat <= c_ord[1]) return(p_ord[1])
  n <- length(c_ord)
  if (stat >= c_ord[n]) return(p_ord[n])
  i <- which(c_ord >= stat)[1]
  if (i <= 1) return(p_ord[1])
  z0 <- qnorm(p_ord[i - 1]); z1 <- qnorm(p_ord[i])
  lam <- (stat - c_ord[i - 1]) / (c_ord[i] - c_ord[i - 1])
  pnorm(z0 + lam * (z1 - z0))
}
