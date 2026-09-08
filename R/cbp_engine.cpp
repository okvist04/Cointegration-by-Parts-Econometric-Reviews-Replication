// cbp_engine.cpp
//
// Rcpp port of fast_engine.jl (the "Cointegration by Parts" Monte Carlo
// engine). The core numerical machinery (fill_prep / window_stat /
// scan_all3) is a direct translation of engine_core.cpp, which was
// validated standalone against a naive brute-force reference to
// floating-point precision (max relative error 1e-13 to 1e-15 across
// all three deterministic cases and N = 1, 3; see verify_engine.R,
// which reruns the same check from R).
//
// RNG: uses R's own RNG (via R::norm_rand() under Rcpp::RNGScope) so
// that R's set.seed() gives full, ordinary reproducibility -- unlike
// Julia's per-replication Xoshiro(seed + r) scheme. Reproducibility
// across replications run via parallel::mclapply() is handled by R's
// own fork-safe RNG stream separation (see R/engine.R).
//
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace Rcpp;

// ---------------------------------------------------------------------
// Deterministic specification (case: 0 = n, 1 = c, 2 = ct)
// ---------------------------------------------------------------------

static inline int det_dim(int case_code) {
    if (case_code == 0) return 0;
    if (case_code == 1) return 1;
    if (case_code == 2) return 2;
    stop("case_code must be 0 (n), 1 (c), or 2 (ct)");
}

static inline int pmax_rule(int n, double mult) {
    int v = (int) std::floor(mult * std::pow((double)n / 100.0, 0.25));
    return std::max(0, v);
}

// [[Rcpp::export]]
int min_window_cpp(int T, double r0, int floor_obs) {
    return std::max((int) std::floor(r0 * T), floor_obs);
}

// [[Rcpp::export]]
int pmax_rule_cpp(int n, double mult = 12.0) {
    return pmax_rule(n, mult);
}

// ---------------------------------------------------------------------
// In-place Cholesky kernels (see engine_core.cpp for derivation notes)
// ---------------------------------------------------------------------

static inline bool chol_upper(std::vector<double>& M, int ld, int n) {
    for (int i = 0; i < n; i++) {
        double s = M[i * ld + i];
        for (int l = 0; l < i; l++) s -= M[l * ld + i] * M[l * ld + i];
        if (!(s > 0)) return false;
        double d = std::sqrt(s);
        M[i * ld + i] = d;
        for (int j = i + 1; j < n; j++) {
            double t = M[i * ld + j];
            for (int l = 0; l < i; l++) t -= M[l * ld + i] * M[l * ld + j];
            M[i * ld + j] = t / d;
        }
    }
    return true;
}

static inline void chol_solve(const std::vector<double>& M, int ld, std::vector<double>& b, int n) {
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int l = 0; l < i; l++) s -= M[l * ld + i] * b[l];
        b[i] = s / M[i * ld + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int l = i + 1; l < n; l++) s -= M[i * ld + l] * b[l];
        b[i] = s / M[i * ld + i];
    }
}

// ---------------------------------------------------------------------
// Prep / Buf (see engine_core.cpp for full derivation notes)
// ---------------------------------------------------------------------

struct Prep {
    int T, N, case_code;
    int k, kk;
    int pmax, nidx, npairs;
    std::vector<double> Szz, Szy, S;
    std::vector<int> pid;
    std::vector<double> Z, Cscratch;

    Prep(int T_, int N_, int case_code_, double pmax_mult)
        : T(T_), N(N_), case_code(case_code_) {
        k = det_dim(case_code) + N;
        kk = k + 1;
        pmax = pmax_rule(T, pmax_mult);
        nidx = pmax + 2;
        npairs = nidx * (nidx + 1) / 2;
        pid.assign(nidx * nidx, 0);
        int c = 0;
        for (int a = 0; a < nidx; a++)
            for (int b = a; b < nidx; b++) { pid[a * nidx + b] = c; pid[b * nidx + a] = c; c++; }
        Szz.assign((size_t)k * k * T, 0.0);
        Szy.assign((size_t)k * T, 0.0);
        S.assign((size_t)kk * kk * npairs * T, 0.0);
        Z.assign((size_t)T * k, 0.0);
        Cscratch.assign((size_t)kk * nidx, 0.0);
    }
    inline double& Szz_at(int i, int j, int t) { return Szz[((size_t)t * k + i) * k + j]; }
    inline double& Szy_at(int i, int t)        { return Szy[(size_t)t * k + i]; }
    inline double& S_at(int i, int j, int p, int t) { return S[(((size_t)t * npairs + p) * kk + i) * kk + j]; }
    inline double& Z_at(int t, int i) { return Z[(size_t)t * k + i]; }
};

static void fill_prep(Prep& P, const std::vector<double>& y, const std::vector<double>& x) {
    int T = P.T, k = P.k, kk = P.kk, N = P.N;
    int kd = det_dim(P.case_code);
    for (int t = 0; t < T; t++) {
        if (kd >= 1) P.Z_at(t, 0) = 1.0;
        if (kd == 2) P.Z_at(t, 1) = (double)(t + 1) / T;
        for (int j = 0; j < N; j++) P.Z_at(t, kd + j) = x[(size_t)t * N + j];
    }
    for (int t = 0; t < T; t++) {
        for (int i = 0; i < k; i++) {
            double zi = P.Z_at(t, i);
            for (int j = 0; j < k; j++) {
                double prev = (t > 0) ? P.Szz_at(i, j, t - 1) : 0.0;
                P.Szz_at(i, j, t) = prev + zi * P.Z_at(t, j);
            }
            double prevy = (t > 0) ? P.Szy_at(i, t - 1) : 0.0;
            P.Szy_at(i, t) = prevy + zi * y[t];
        }
    }
    for (int t = 0; t < T; t++) {
        std::fill(P.Cscratch.begin(), P.Cscratch.end(), 0.0);
        if (t >= 1) {
            P.Cscratch[0 * P.nidx + 0] = y[t - 1];
            for (int i = 0; i < k; i++) P.Cscratch[(1 + i) * P.nidx + 0] = P.Z_at(t - 1, i);
        }
        for (int j = 0; j <= P.pmax; j++) {
            int tt = t - j;
            if (tt >= 1) {
                int col = 1 + j;
                P.Cscratch[0 * P.nidx + col] = y[tt] - y[tt - 1];
                for (int i = 0; i < k; i++)
                    P.Cscratch[(1 + i) * P.nidx + col] = P.Z_at(tt, i) - P.Z_at(tt - 1, i);
            }
        }
        for (int a = 0; a < P.nidx; a++) {
            for (int b = a; b < P.nidx; b++) {
                int p = P.pid[a * P.nidx + b];
                for (int i = 0; i < kk; i++) {
                    double cia = P.Cscratch[i * P.nidx + a];
                    for (int j = 0; j < kk; j++) {
                        double prev = (t > 0) ? P.S_at(i, j, p, t - 1) : 0.0;
                        P.S_at(i, j, p, t) = prev + cia * P.Cscratch[j * P.nidx + b];
                    }
                }
            }
        }
    }
}

struct Buf {
    int k, kk, pmax, Gld;
    std::vector<double> A, b, a, G, bet, w;
    explicit Buf(const Prep& P) : k(P.k), kk(P.kk), pmax(P.pmax) {
        A.assign((size_t)k * k, 0.0); b.assign(k, 0.0); a.assign(kk, 0.0);
        Gld = pmax + 2; G.assign((size_t)Gld * Gld, 0.0);
        bet.assign(pmax + 1, 0.0); w.assign(pmax + 1, 0.0);
    }
};

struct WinResult { double tstat; int lag; };

static WinResult window_stat(Prep& P, int s, int e, Buf& B, double pmax_mult, double dof_ratio, bool aic) {
    int n = e - s + 1;
    int k = P.k, kk = P.kk;
    for (int i = 0; i < k; i++) {
        double prevb = (s > 0) ? P.Szy_at(i, s - 1) : 0.0;
        B.b[i] = P.Szy_at(i, e) - prevb;
        for (int j = 0; j < k; j++) {
            double prevA = (s > 0) ? P.Szz_at(i, j, s - 1) : 0.0;
            B.A[i * k + j] = P.Szz_at(i, j, e) - prevA;
        }
    }
    if (!chol_upper(B.A, k, k)) return {NA_REAL, -1};
    chol_solve(B.A, k, B.b, k);
    B.a[0] = 1.0;
    for (int i = 0; i < k; i++) B.a[1 + i] = -B.b[i];

    int pw = std::min(P.pmax, pmax_rule(n, pmax_mult));
    while (pw > 0 && (n - 1 - pw) < std::max(pw + 3, (int)std::ceil(dof_ratio * (pw + 1)))) pw -= 1;
    int neff = n - 1 - pw;
    if (neff < 3) return {NA_REAL, -1};
    int m = pw + 1, dim = m + 1, lo = s + pw;

    std::vector<double>& G = B.G;
    int Gld = B.Gld;
    for (int u = 0; u < dim; u++) {
        int ia = (u == 0) ? 0 : ((u == dim - 1) ? 1 : u + 1);
        for (int v = u; v < dim; v++) {
            int ib = (v == 0) ? 0 : ((v == dim - 1) ? 1 : v + 1);
            int pp = P.pid[ia * P.nidx + ib];
            double acc = 0.0;
            for (int i = 0; i < kk; i++) {
                double ai = B.a[i];
                if (ai == 0.0) continue;
                for (int j = 0; j < kk; j++) {
                    double hi = P.S_at(i, j, pp, e);
                    double lov = (lo >= 0) ? P.S_at(i, j, pp, lo) : 0.0;
                    acc += ai * B.a[j] * (hi - lov);
                }
            }
            G[u * Gld + v] = acc; G[v * Gld + u] = acc;
        }
    }
    double yy = G[(dim - 1) * Gld + (dim - 1)];
    if (!(yy > 0)) return {NA_REAL, -1};
    if (!chol_upper(G, Gld, dim)) return {NA_REAL, -1};

    double best_ic = std::numeric_limits<double>::infinity();
    int best_p = -1; double best_rss = NA_REAL, acc = 0.0;
    for (int p = 0; p <= pw; p++) {
        double gpd = G[p * Gld + (dim - 1)];
        acc += gpd * gpd;
        double rss = yy - acc;
        int kp = p + 1;
        if (!(rss > 0 && neff > kp)) continue;
        double pen = aic ? 2.0 : std::log((double)neff);
        double ic = neff * std::log(rss / neff) + pen * kp;
        if (ic < best_ic) { best_ic = ic; best_p = p; best_rss = rss; }
    }
    if (best_p < 0) return {NA_REAL, -1};

    int j = best_p + 1;
    for (int i = j - 1; i >= 0; i--) {
        double ssum = G[i * Gld + (dim - 1)];
        for (int l = i + 1; l < j; l++) ssum -= G[i * Gld + l] * B.bet[l];
        B.bet[i] = ssum / G[i * Gld + i];
    }
    double gam = B.bet[0];
    B.w[0] = 1.0 / G[0 * Gld + 0];
    for (int i = 1; i < j; i++) {
        double ssum = 0.0;
        for (int l = 0; l < i; l++) ssum -= G[l * Gld + i] * B.w[l];
        B.w[i] = ssum / G[i * Gld + i];
    }
    double v11 = 0.0; for (int i = 0; i < j; i++) v11 += B.w[i] * B.w[i];
    double sig2 = best_rss / (neff - j);
    if (!(sig2 > 0 && v11 > 0)) return {NA_REAL, -1};
    double se = std::sqrt(sig2 * v11);
    if (!(se > 0)) return {NA_REAL, -1};
    return { gam / se, best_p };
}

struct Scan3Result { double gieg, fieg, bieg; int gs, ge, fe, bs; };

static Scan3Result scan_all3(Prep& P, int min_obs, Buf& B, double pmax_mult, double dof_ratio, bool aic) {
    int T = P.T;
    double g = std::numeric_limits<double>::infinity(); int gs = -1, ge = -1;
    double f = std::numeric_limits<double>::infinity(); int fe = -1;
    double bk = std::numeric_limits<double>::infinity(); int bs = -1;
    for (int s = 0; s <= T - min_obs; s++) {
        for (int e = s + min_obs - 1; e < T; e++) {
            WinResult r = window_stat(P, s, e, B, pmax_mult, dof_ratio, aic);
            if (!R_finite(r.tstat)) continue;
            if (r.tstat < g) { g = r.tstat; gs = s; ge = e; }
            if (s == 0 && r.tstat < f) { f = r.tstat; fe = e; }
            if (e == T - 1 && r.tstat < bk) { bk = r.tstat; bs = s; }
        }
    }
    return { R_finite(g) ? g : NA_REAL, R_finite(f) ? f : NA_REAL, R_finite(bk) ? bk : NA_REAL, gs, ge, fe, bs };
}

// ---------------------------------------------------------------------
// Data generating processes (R's RNG, so set.seed() controls them)
// ---------------------------------------------------------------------

static void dgp_null(std::vector<double>& y, std::vector<double>& x, int T, int N,
                      double drift_x, double drift_u) {
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int t = 0; t < T; t++) { acc += R::norm_rand() + drift_x; x[(size_t)t * N + j] = acc; }
    }
    double acc = 0.0;
    for (int t = 0; t < T; t++) { acc += R::norm_rand() + drift_u; y[t] = acc; }
}

// direction_code: 0 = forward (I(0) then I(1)), 1 = reverse (I(1) then I(0))
static void dgp_break(std::vector<double>& y, std::vector<double>& x, int T, int N,
                       double tau0, int direction_code, double rho) {
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int t = 0; t < T; t++) { acc += R::norm_rand(); x[(size_t)t * N + j] = acc; }
    }
    int tb = (int) std::floor(tau0 * T);
    double e = 0.0;
    for (int t = 0; t < T; t++) {
        bool stationary = (direction_code == 0) ? (t < tb) : (t >= tb); // 0-based: t<tb <=> (t+1)<=tb in 1-based
        e = stationary ? rho * e + R::norm_rand() : e + R::norm_rand();
        double s = 0.0; for (int j = 0; j < N; j++) s += x[(size_t)t * N + j];
        y[t] = s + e;
    }
}

static void dgp_window(std::vector<double>& y, std::vector<double>& x, int T, int N,
                        double tau1, double tau2, double rho, int& t1_out, int& t2_out) {
    for (int j = 0; j < N; j++) {
        double acc = 0.0;
        for (int t = 0; t < T; t++) { acc += R::norm_rand(); x[(size_t)t * N + j] = acc; }
    }
    int t1 = (int) std::floor(tau1 * T);
    int t2 = (int) std::floor(tau2 * T);
    t1_out = t1; t2_out = t2;
    double e = 0.0;
    for (int t = 0; t < T; t++) {
        bool stationary = (t >= t1) && (t < t2); // 0-based: (t>t1)&&(t<=t2) in 1-based -> t+1>t1 && t+1<=t2 -> t>=t1 && t<t2
        e = stationary ? rho * e + R::norm_rand() : e + R::norm_rand();
        double s = 0.0; for (int j = 0; j < N; j++) s += x[(size_t)t * N + j];
        y[t] = s + e;
    }
}

// ---------------------------------------------------------------------
// Naive reference (for verify_engine.R); small dense OLS via Cholesky
// ---------------------------------------------------------------------

static bool ols_solve(const std::vector<double>& Zmat, int n, int p, const std::vector<double>& yv,
                       std::vector<double>& beta_out, std::vector<double>& resid_out, double& v11_out) {
    std::vector<double> XtX((size_t)p * p, 0.0), Xty(p, 0.0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int t = 0; t < n; t++) s += Zmat[(size_t)t * p + i] * Zmat[(size_t)t * p + j];
            XtX[i * p + j] = s;
        }
        double sy = 0.0; for (int t = 0; t < n; t++) sy += Zmat[(size_t)t * p + i] * yv[t];
        Xty[i] = sy;
    }
    std::vector<double> M = XtX;
    if (!chol_upper(M, p, p)) return false;
    std::vector<double> b = Xty;
    chol_solve(M, p, b, p);
    beta_out = b;
    resid_out.assign(n, 0.0);
    for (int t = 0; t < n; t++) {
        double pred = 0.0; for (int i = 0; i < p; i++) pred += Zmat[(size_t)t * p + i] * beta_out[i];
        resid_out[t] = yv[t] - pred;
    }
    std::vector<double> w(p, 0.0); w[0] = 1.0;
    chol_solve(M, p, w, p);
    v11_out = w[0];
    return true;
}

struct AdfResult { double tstat; int lag; };

static AdfResult adf_naive(const std::vector<double>& series, double pmax_mult, bool common_sample, bool aic) {
    int n0 = (int) series.size();
    if (n0 < 10) return {NA_REAL, -1};
    int pm = pmax_rule(n0, pmax_mult);
    std::vector<double> dy(n0 - 1), ylag(n0 - 1);
    for (int i = 0; i < n0 - 1; i++) { dy[i] = series[i + 1] - series[i]; ylag[i] = series[i]; }
    if (common_sample) while (pm > 0 && ((int)dy.size() - pm) < (pm + 3)) pm -= 1;
    double best_ic = std::numeric_limits<double>::infinity();
    double best_t = NA_REAL; int best_p = -1;
    for (int p = 0; p <= pm; p++) {
        int start = common_sample ? pm : p;
        int n = (int)dy.size() - start;
        if (n <= 0) continue;
        int cols = p + 1;
        if (n <= cols) continue;
        std::vector<double> X((size_t)n * cols, 0.0), Y(n, 0.0);
        for (int t = 0; t < n; t++) {
            Y[t] = dy[start + t];
            X[(size_t)t * cols + 0] = ylag[start + t];
            for (int jl = 1; jl <= p; jl++) X[(size_t)t * cols + jl] = dy[start + t - jl];
        }
        std::vector<double> beta, resid; double v11;
        if (!ols_solve(X, n, cols, Y, beta, resid, v11)) continue;
        double rss = 0.0; for (double rv : resid) rss += rv * rv;
        if (!(rss > 0)) continue;
        double sig2 = rss / (n - cols);
        if (!(sig2 > 0 && v11 > 0)) continue;
        double tst = beta[0] / std::sqrt(sig2 * v11);
        double ic = n * std::log(rss / n) + (aic ? 2.0 : std::log((double)n)) * cols;
        if (ic < best_ic) { best_ic = ic; best_t = tst; best_p = p; }
    }
    return { best_t, best_p };
}

static AdfResult window_stat_naive(const std::vector<double>& y, const std::vector<double>& x,
                                    int s, int e, int case_code, int N,
                                    double pmax_mult, bool common_sample) {
    int n = e - s + 1;
    int kd = det_dim(case_code);
    int p = kd + N;
    std::vector<double> Z((size_t)n * p, 0.0), yv(n, 0.0);
    for (int t = 0; t < n; t++) {
        int col = 0;
        if (kd >= 1) Z[(size_t)t * p + col++] = 1.0;
        if (kd == 2) Z[(size_t)t * p + col++] = (double)(t + 1) / n;
        for (int j = 0; j < N; j++) Z[(size_t)t * p + col + j] = x[(size_t)(s + t) * N + j];
        yv[t] = y[s + t];
    }
    std::vector<double> beta, resid; double v11;
    if (!ols_solve(Z, n, p, yv, beta, resid, v11)) return {NA_REAL, -1};
    return adf_naive(resid, pmax_mult, common_sample, false);
}

// ---------------------------------------------------------------------
// Gregory-Hansen benchmark (model C: level shift), naive but adequate
// ---------------------------------------------------------------------

// [[Rcpp::export]]
double gh_stat_cpp(NumericVector y_r, NumericMatrix x_r, double trim = 0.15, double pmax_mult = 12.0) {
    int T = y_r.size(), N = x_r.ncol();
    std::vector<double> y(y_r.begin(), y_r.end());
    // Julia (1-based): lo = max(2, floor(trim*T)), hi = min(T-1, ceil((1-trim)*T)).
    // 0-based b_cpp = b_julia - 1, with the dummy condition t_cpp > b_cpp
    // reproducing t_julia > b_julia exactly (checked against a worked example).
    int lo = std::max(1, (int)std::floor(trim * T) - 1);   // 0-based
    int hi = std::min(T - 2, (int)std::ceil((1 - trim) * T) - 1);
    int stride = std::max(1, T / 150);
    double best = std::numeric_limits<double>::infinity();
    int p = 2 + N;
    std::vector<double> Z((size_t)T * p, 0.0);
    for (int t = 0; t < T; t++) {
        Z[(size_t)t * p + 0] = 1.0;
        for (int j = 0; j < N; j++) Z[(size_t)t * p + 2 + j] = x_r(t, j);
    }
    for (int b = lo; b <= hi; b += stride) {
        for (int t = 0; t < T; t++) Z[(size_t)t * p + 1] = (t > b) ? 1.0 : 0.0;
        std::vector<double> beta, resid; double v11;
        if (!ols_solve(Z, T, p, y, beta, resid, v11)) continue;
        AdfResult r = adf_naive(resid, pmax_mult, true, false);
        if (R_finite(r.tstat) && r.tstat < best) best = r.tstat;
    }
    return R_finite(best) ? best : NA_REAL;
}

// ---------------------------------------------------------------------
// Per-replication wrappers: DGP + fill_prep + scan_all3 (+ static EG)
// in one call, so R only pays one function-call boundary per
// replication when looping via parallel::mclapply().
// ---------------------------------------------------------------------

// [[Rcpp::export]]
List rep_null_cpp(int T, int N, int case_code, int min_obs,
                   double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false,
                   double drift = 0.0) {
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    dgp_null(y, x, T, N, drift, drift);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1 // 1-based for R
    );
}

// [[Rcpp::export]]
List rep_break_cpp(int T, int N, int case_code, int min_obs,
                    double tau0, int direction_code, double rho,
                    double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    dgp_break(y, x, T, N, tau0, direction_code, rho);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1
    );
}

// [[Rcpp::export]]
List rep_window_cpp(int T, int N, int case_code, int min_obs,
                     double tau1, double tau2, double rho,
                     double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    int t1, t2;
    dgp_window(y, x, T, N, tau1, tau2, rho, t1, t2);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    // t1, t2 are floor(tau*T) COUNTS (not positions), used identically in
    // Julia and R's scoring formula (inside = gs >= t1+1 & ge <= t2) --
    // do NOT shift these to 1-based, only gs/ge/fe/bs (actual indices) are.
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1,
        _["t1"] = t1, _["t2"] = t2
    );
}

// [[Rcpp::export]]
List rep_null_gh_cpp(int T, int N, int case_code, int min_obs,
                      double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    // like rep_null_cpp but also computes the Gregory-Hansen statistic on
    // the same draw, for run_size_power.R's GH benchmark
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    dgp_null(y, x, T, N, 0.0, 0.0);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    NumericMatrix xm(T, N);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) xm(t, j) = x[(size_t)t * N + j];
    NumericVector yv(y.begin(), y.end());
    double gh = gh_stat_cpp(yv, xm, 0.15, pmax_mult);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat, _["gh"] = gh,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1
    );
}

// [[Rcpp::export]]
List rep_break_gh_cpp(int T, int N, int case_code, int min_obs,
                       double tau0, int direction_code, double rho,
                       double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    dgp_break(y, x, T, N, tau0, direction_code, rho);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    NumericMatrix xm(T, N);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) xm(t, j) = x[(size_t)t * N + j];
    NumericVector yv(y.begin(), y.end());
    double gh = gh_stat_cpp(yv, xm, 0.15, pmax_mult);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat, _["gh"] = gh,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1
    );
}

// [[Rcpp::export]]
List rep_window_gh_cpp(int T, int N, int case_code, int min_obs,
                        double tau1, double tau2, double rho,
                        double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    // like rep_window_cpp but also computes the Gregory-Hansen statistic
    // on the same draw, since run_size_power.jl applies the GH benchmark
    // uniformly to every power configuration, window included.
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    int t1, t2;
    dgp_window(y, x, T, N, tau1, tau2, rho, t1, t2);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    WinResult eg = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    NumericMatrix xm(T, N);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) xm(t, j) = x[(size_t)t * N + j];
    NumericVector yv(y.begin(), y.end());
    double gh = gh_stat_cpp(yv, xm, 0.15, pmax_mult);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg, _["eg"] = eg.tstat, _["gh"] = gh,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1,
        _["t1"] = t1, _["t2"] = t2
    );
}

// ---------------------------------------------------------------------
// Direct scan on user-supplied data (application, verification, and
// the recursive/second-stage subsample logic)
// ---------------------------------------------------------------------

// [[Rcpp::export]]
List scan_all3_data_cpp(NumericVector y_r, NumericMatrix x_r, int case_code, int min_obs,
                         double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    int T = y_r.size(), N = x_r.ncol();
    std::vector<double> y(y_r.begin(), y_r.end());
    std::vector<double> x((size_t)T * N, 0.0);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) x[(size_t)t * N + j] = x_r(t, j);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    Scan3Result r = scan_all3(P, min_obs, B, pmax_mult, dof_ratio, aic);
    return List::create(
        _["fieg"] = r.fieg, _["bieg"] = r.bieg, _["gieg"] = r.gieg,
        _["fe"] = r.fe + 1, _["bs"] = r.bs + 1, _["gs"] = r.gs + 1, _["ge"] = r.ge + 1
    );
}

// [[Rcpp::export]]
double eg_full_data_cpp(NumericVector y_r, NumericMatrix x_r, int case_code,
                         double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    int T = y_r.size(), N = x_r.ncol();
    std::vector<double> y(y_r.begin(), y_r.end());
    std::vector<double> x((size_t)T * N, 0.0);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) x[(size_t)t * N + j] = x_r(t, j);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    WinResult r = window_stat(P, 0, T - 1, B, pmax_mult, dof_ratio, aic);
    return r.tstat;
}

// [[Rcpp::export]]
List adf_series_cpp(NumericVector v_r, bool trend, double pmax_mult = 12.0) {
    // ADF pre-test on a single series with a constant (and optionally a
    // trend) removed first -- the integration-order pre-test.
    int n = v_r.size();
    std::vector<double> Z((size_t)n * (trend ? 2 : 1), 0.0);
    int p = trend ? 2 : 1;
    for (int t = 0; t < n; t++) {
        Z[(size_t)t * p + 0] = 1.0;
        if (trend) Z[(size_t)t * p + 1] = (double)(t + 1) / n;
    }
    std::vector<double> v(v_r.begin(), v_r.end());
    std::vector<double> beta, resid; double v11;
    if (!ols_solve(Z, n, p, v, beta, resid, v11)) return List::create(_["tstat"] = NA_REAL, _["lag"] = -1);
    AdfResult r = adf_naive(resid, pmax_mult, true, false);
    return List::create(_["tstat"] = r.tstat, _["lag"] = r.lag);
}

// ---------------------------------------------------------------------
// Verification-only exports: single window, fast vs naive, brute-force
// scan -- used by verify_engine.R, kept slow-but-simple deliberately.
// ---------------------------------------------------------------------

// [[Rcpp::export]]
List window_stat_cpp(NumericVector y_r, NumericMatrix x_r, int case_code, int s1, int e1,
                      double pmax_mult = 12.0, double dof_ratio = 0.0, bool aic = false) {
    // s1, e1 are 1-based (R convention) inclusive window bounds
    int T = y_r.size(), N = x_r.ncol();
    std::vector<double> y(y_r.begin(), y_r.end());
    std::vector<double> x((size_t)T * N, 0.0);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) x[(size_t)t * N + j] = x_r(t, j);
    Prep P(T, N, case_code, pmax_mult);
    fill_prep(P, y, x);
    Buf B(P);
    WinResult r = window_stat(P, s1 - 1, e1 - 1, B, pmax_mult, dof_ratio, aic);
    return List::create(_["tstat"] = r.tstat, _["lag"] = r.lag);
}

// [[Rcpp::export]]
List window_stat_naive_cpp(NumericVector y_r, NumericMatrix x_r, int case_code, int s1, int e1,
                            double pmax_mult = 12.0, bool common_sample = true) {
    int T = y_r.size(), N = x_r.ncol();
    std::vector<double> y(y_r.begin(), y_r.end());
    std::vector<double> x((size_t)T * N, 0.0);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) x[(size_t)t * N + j] = x_r(t, j);
    AdfResult r = window_stat_naive(y, x, s1 - 1, e1 - 1, case_code, N, pmax_mult, common_sample);
    return List::create(_["tstat"] = r.tstat, _["lag"] = r.lag);
}

// [[Rcpp::export]]
List dgp_null_cpp(int T, int N, double drift = 0.0) {
    RNGScope scope;
    std::vector<double> y(T, 0.0), x((size_t)T * N, 0.0);
    dgp_null(y, x, T, N, drift, drift);
    NumericMatrix xm(T, N);
    for (int t = 0; t < T; t++) for (int j = 0; j < N; j++) xm(t, j) = x[(size_t)t * N + j];
    return List::create(_["y"] = NumericVector(y.begin(), y.end()), _["x"] = xm);
}
