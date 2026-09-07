// adf_scan_fast.cpp
//
// Fast C++ replacement for the hot loops in generate_cv.R:
//   residual_df_window(), adf_noconst_fast(),
//   scan_expanding_windows(), scan_all_subsamples()
//
// Statistical logic is kept IDENTICAL to the R version (same lag grid,
// same BIC/AIC selection, same intercept handling). Two things differ,
// both purely mechanical / non-statistical:
//
//   1. Per-window result tables (start/stop/tstat/lag for every window)
//      are never materialized -- the null-calibration pipeline only
//      ever uses the running minimum t-stat, so we track that directly.
//      This alone removes a huge amount of allocation for "all", which
//      has ~T^2/2 windows.
//
//   2. Each candidate lag order is fit with ONE Cholesky factorization
//      of X'X plus two triangular solves, instead of the R version's
//      separate qr()-rank-check + lm.fit() (a second QR internally) +
//      a full matrix inverse when only one element of it is used.
//
// Two exported entry points:
//   - minstat_scan_cpp(): returns only the minimum statistic. Used by
//     the Monte Carlo calibration (generate_cv.R), where only the value
//     is needed across 10,000s of replications and per-window detail
//     would be wasted allocation.
//   - scan_full_cpp(): also returns the location of the extremizing
//     window and, optionally, the full scan path. Used for one-off
//     application to real data (run_scan_test.R), where that detail is
//     exactly what's wanted and the cost of keeping it is negligible.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

struct AdfResult {
  double tstat;
  int lag;
};

// Fit y = X*beta via a single Cholesky factorization of X'X.
// Returns false if X'X is not positive definite (rank-deficient /
// too few effective observations for the given lag order) -- this is
// the direct analogue of the "qr(X)$rank < ncol(X)" guard in the R code.
static bool chol_fit(const arma::mat &X, const arma::vec &y,
                      arma::vec &beta, arma::vec &resid,
                      double &inv00) {
  arma::mat XtX = X.t() * X;
  arma::mat R;
  if (!arma::chol(R, XtX)) {
    return false;
  }
  arma::vec Xty = X.t() * y;

  // Solve R^T R beta = Xty via two triangular solves (avoids computing
  // the full inverse just to get beta).
  arma::vec w = arma::solve(arma::trimatl(R.t()), Xty);
  beta = arma::solve(arma::trimatu(R), w);
  resid = y - X * beta;

  // (X'X)^{-1}[0,0], needed for se_phi, via one more triangular solve
  // instead of a full k x k matrix inverse.
  arma::vec e1(X.n_cols, arma::fill::zeros);
  e1(0) = 1.0;
  arma::vec x = arma::solve(arma::trimatl(R.t()), e1);
  inv00 = arma::dot(x, x);
  return true;
}

// Equivalent of adf_noconst_fast(): BIC/AIC lag-order selection on a
// series with no deterministic terms.
static AdfResult adf_noconst_fast_cpp(const arma::vec &series,
                                       int max_lag_input,
                                       const std::string &criterion) {
  int T = series.n_elem;
  AdfResult best;
  best.tstat = NA_REAL;
  best.lag = -1;

  if (T < 10) return best;

  int max_lag = max_lag_input;
  if (max_lag < 0) {
    // Same default rule as the R version, applied to THIS window's
    // own length -- matches the original per-window dynamic behavior.
    max_lag = (int) std::floor(12.0 * std::pow(T / 100.0, 0.25));
  }
  if (max_lag < 0) max_lag = 0;

  int nd = T - 1; // length(dy)
  arma::vec dy = series.subvec(1, T - 1) - series.subvec(0, T - 2);
  arma::vec ylag = series.subvec(0, T - 2);

  double best_ic = R_PosInf;

  for (int p = 0; p <= max_lag; p++) {
    int n = nd - p;
    if (n <= 0) continue;

    arma::vec Y = dy.subvec(p, nd - 1);
    arma::vec yreg = ylag.subvec(p, nd - 1);

    arma::mat X;
    if (p == 0) {
      X = yreg;
    } else {
      arma::mat Xlags(n, p);
      for (int j = 1; j <= p; j++) {
        Xlags.col(j - 1) = dy.subvec(p - j, nd - 1 - j);
      }
      X = arma::join_horiz(yreg, Xlags);
    }

    int k = X.n_cols;
    if (n <= k) continue;

    arma::vec beta, resid;
    double inv00;
    if (!chol_fit(X, Y, beta, resid, inv00)) continue;

    double rss = arma::dot(resid, resid);
    if (rss <= 0) continue;

    double sigma2 = rss / (n - k);
    if (inv00 <= 0) continue;
    double se_phi = std::sqrt(sigma2 * inv00);
    if (se_phi <= 0) continue;

    double tstat = beta(0) / se_phi;

    double ic;
    if (criterion == "aic") {
      ic = n * std::log(rss / n) + 2.0 * k;
    } else {
      ic = n * std::log(rss / n) + std::log((double) n) * k;
    }

    if (ic < best_ic) {
      best_ic = ic;
      best.tstat = tstat;
      best.lag = p;
    }
  }

  return best;
}

// Equivalent of residual_df_window(): first-stage OLS (with optional
// intercept), then ADF-without-constant on the residuals.
static AdfResult residual_df_window_cpp(const arma::vec &y_sub,
                                         const arma::mat &X_sub,
                                         bool include_intercept,
                                         int max_lag) {
  AdfResult none;
  none.tstat = NA_REAL;
  none.lag = -1;

  int n = y_sub.n_elem;
  if (n < 10) return none;

  arma::mat X;
  if (include_intercept) {
    arma::vec ones(n, arma::fill::ones);
    X = arma::join_horiz(ones, X_sub);
  } else {
    X = X_sub;
  }

  arma::vec beta, resid;
  double inv00;
  if (!chol_fit(X, y_sub, beta, resid, inv00)) return none;

  return adf_noconst_fast_cpp(resid, max_lag, "bic");
}

// [[Rcpp::export]]
double minstat_scan_cpp(NumericVector y_r, NumericMatrix X_r,
                         std::string scan_type,
                         int min_obs, int min_start, int min_window_obs,
                         bool include_intercept, int max_lag) {

  arma::vec y(y_r.begin(), y_r.size(), false);
  arma::mat X(X_r.begin(), X_r.nrow(), X_r.ncol(), false);
  int T = y.n_elem;

  double min_t = R_PosInf;
  bool found = false;

  if (scan_type == "fw" || scan_type == "bw") {

    int n_windows = T - min_obs + 1;

    for (int i = 1; i <= n_windows; i++) {
      int s, e;
      if (scan_type == "fw") {
        s = 1; e = min_obs + i - 1;
      } else {
        s = T - (min_obs + i - 1) + 1; e = T;
      }
      arma::vec y_sub = y.subvec(s - 1, e - 1);
      arma::mat X_sub = X.rows(s - 1, e - 1);
      AdfResult r = residual_df_window_cpp(y_sub, X_sub, include_intercept, max_lag);
      if (std::isfinite(r.tstat) && r.tstat < min_t) { min_t = r.tstat; found = true; }
    }

  } else { // "all"

    int s_max = T - min_window_obs + 1;

    for (int s = min_start; s <= s_max; s++) {
      int e_min = s + min_window_obs - 1;
      for (int e = e_min; e <= T; e++) {
        arma::vec y_sub = y.subvec(s - 1, e - 1);
        arma::mat X_sub = X.rows(s - 1, e - 1);
        AdfResult r = residual_df_window_cpp(y_sub, X_sub, include_intercept, max_lag);
        if (std::isfinite(r.tstat) && r.tstat < min_t) { min_t = r.tstat; found = true; }
      }
    }
  }

  return found ? min_t : NA_REAL;
}

// [[Rcpp::export]]
List scan_full_cpp(NumericVector y_r, NumericMatrix X_r,
                    std::string scan_type,
                    int min_obs, int min_start, int min_window_obs,
                    bool include_intercept, int max_lag,
                    bool return_full_scan) {

  arma::vec y(y_r.begin(), y_r.size(), false);
  arma::mat X(X_r.begin(), X_r.nrow(), X_r.ncol(), false);
  int T = y.n_elem;

  double min_t = R_PosInf;
  bool found = false;
  int best_s = NA_INTEGER, best_e = NA_INTEGER, best_lag = NA_INTEGER;

  // Only used when return_full_scan is true. Reserved up front where the
  // final size is known cheaply; for "all" it can be large (~T^2/2), but
  // this function is meant for one-off application to a single real
  // dataset, not repeated Monte Carlo calls, so that's fine.
  std::vector<int> out_start, out_end, out_lag;
  std::vector<double> out_tstat;

  auto evaluate_window = [&](int s, int e) {
    arma::vec y_sub = y.subvec(s - 1, e - 1);
    arma::mat X_sub = X.rows(s - 1, e - 1);
    AdfResult r = residual_df_window_cpp(y_sub, X_sub, include_intercept, max_lag);

    if (return_full_scan) {
      out_start.push_back(s);
      out_end.push_back(e);
      out_tstat.push_back(r.tstat);
      out_lag.push_back(r.lag);
    }

    if (std::isfinite(r.tstat) && r.tstat < min_t) {
      min_t = r.tstat;
      found = true;
      best_s = s;
      best_e = e;
      best_lag = r.lag;
    }
  };

  if (scan_type == "fw" || scan_type == "bw") {

    int n_windows = T - min_obs + 1;

    for (int i = 1; i <= n_windows; i++) {
      int s, e;
      if (scan_type == "fw") {
        s = 1; e = min_obs + i - 1;
      } else {
        s = T - (min_obs + i - 1) + 1; e = T;
      }
      evaluate_window(s, e);
    }

  } else { // "all"

    int s_max = T - min_window_obs + 1;

    for (int s = min_start; s <= s_max; s++) {
      int e_min = s + min_window_obs - 1;
      for (int e = e_min; e <= T; e++) {
        evaluate_window(s, e);
      }
    }
  }

  List out = List::create(
    Named("min_t") = found ? min_t : NA_REAL,
    Named("best_start") = best_s,
    Named("best_end") = best_e,
    Named("best_lag") = best_lag
  );

  if (return_full_scan) {
    out["start"] = wrap(out_start);
    out["end"] = wrap(out_end);
    out["tstat"] = wrap(out_tstat);
    out["lag"] = wrap(out_lag);
  }

  return out;
}
