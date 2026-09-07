# Cointegration by Parts

Replication package for "Cointegration by Parts" (submitted to *Econometric
Reviews*, Special Issue in Honor of James G. MacKinnon).

This paper develops the **Cointegration by Parts** test — a recursive,
window-based procedure (inspired by recursive unit-root tests used to detect
financial bubbles) for detecting whether a cointegrating relationship breaks
down or is unstable over time. Finite-sample critical values are obtained via
a MacKinnon-style response-surface regression. The method is validated by
simulation and applied to global mean sea level (GMSL) and global mean
land-ocean temperature anomalies (GMTA).

## Repository structure

```
.
├── R/                    Core methodology: simulation, calibration,
│                         response-surface fitting, and the scan test itself.
├── data/                 Raw input data for the empirical application.
├── analysis/             The empirical application script (GMSL vs. GMTA).
├── output/               All generated figures and calibration outputs.
└── CointByParts.Rproj    RStudio project file.
```

### `R/` — core methodology

| File | Purpose |
|---|---|
| `adf_scan_fast.cpp` | Rcpp/Armadillo engine for the recursive ADF scan (`scan_full_cpp`, used by the real-data test; `minstat_scan_cpp`, used by the null-calibration simulation, which only needs the running minimum). |
| `generate_cv.R` | Simulation and calibration machinery: data-generating processes under the null and under a break (`simulate_pair_null`, `simulate_window_break`, `simulate_pair_break`), the fast (C++) and pure-R scan paths, and `calibrate_null_critical_values_parallel()`, the top-level function that builds the null distribution of the test statistic across a grid of `(T, n_regressors, scan_type)`. |
| `run_generation_cv.R` | Runs the full calibration grid (`T` from 50 to 1000, 1–6 regressors, all three scan types, 10,000 replications each) and writes `output/cv_output/critical_values_null_full_grid.csv`. Checkpointed per combination — safe to interrupt and re-run. |
| `validate_fast_vs_r.R` | **Run this before trusting the calibration grid.** Compares the fast C++ path against the pure-R reference implementation on identical simulated data; writes `output/cv_output/validation_fast_vs_r.csv`. Differences should be at floating-point tolerance (< 1e-6). |
| `fit_response_surface.R` | Fits the response surface (`crit_val ~ b0 + b1/T + b2/T² [+ b3/T³]`) to the calibration grid, per `(scan_type, n_regressors, prob)`, and provides `get_critical_value()` for looking up a critical value at any `T`. Writes `output/cv_output/response_surface_coefficients.csv`. |
| `run_scan_test.R` | The test itself: `run_scan_test()` / `test_cointegration()` run FIEG/BIEG/GIEG on real data and return a `scan_test` object; `add_critical_values()` attaches significance from the fitted response surface; `print.scan_test()` / `plot.scan_test()` are the display methods. |
| `motivate_simulated_break.R` | Motivating simulation example: generates data with a coint → non-coint → coint regime pattern and shows that GIEG(r0) correctly isolates the cointegrated segment. Produces `output/figures/sim_breaks.eps` and `output/figures/t_stat_sim_breaks.eps`. |

### `data/`

- `CSIRO_Recons_gmsl_yr_2019.csv` — CSIRO reconstructed global mean sea level.
- `HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv` — HadCRUT5 global annual temperature anomaly (native 1961–1990 baseline; rebased to 1850–1900 in the analysis script).

### `analysis/empirical_application.R`

The full empirical pipeline, in order:

1. Loads and aligns GMSL and GMTA to a common annual sample; rebases HadCRUT5 to a 1850–1900 baseline.
2. Plots the two series on a dual axis → `output/figures/GMST_GMTA.eps`.
3. Univariate ADF tests (drift, BIC-selected lags) confirming both series are non-stationary in levels. Almost a replication of Schmith, Johansen & Thejll (2012)'s ADF specification (constant only, 2 fixed lags), but without 2 fixed lags.
4. Runs FIEG(r0), BIEG(r0,1), and GIEG(r0) on the full sample, with critical values from the fitted response surface.
5. A **recursive break search**: after finding a significant GIEG window, the sample is split into the two remaining subsamples (excluding that window) and re-tested, recursing until a subsample falls below the minimum usable size.
6. Extracts all significant windows found across the recursive search and plots them shaded on the original series → `output/figures/cointegration_windows.eps`.

## Mapping to the paper

| Paper item | Produced by |
|---|---|
| Response-surface coefficients table | `R/fit_response_surface.R` → `output/cv_output/response_surface_coefficients.csv` |
| Motivating simulation figure(s) | `R/motivate_simulated_break.R` → `output/figures/sim_breaks.eps`, `output/figures/t_stat_sim_breaks.eps` |
| GMSL / GMTA series figure | `analysis/empirical_application.R` (Step 2) → `output/figures/GMST_GMTA.eps` |
| Univariate ADF pre-tests (order of integration) | `analysis/empirical_application.R` (Step 3) |
| Full-sample FIEG / BIEG / GIEG results table | `analysis/empirical_application.R` (Step 4) |
| Subsample (recursive break search) results table(s) | `analysis/empirical_application.R` (Step 5) |
| Cointegrated-windows figure | `analysis/empirical_application.R` (Step 6) → `output/figures/cointegration_windows.eps` |

*(Table/figure numbers are left for you to fill in against the current paper draft — the mapping above is by content, not by number, since those may still shift before submission.)*

## Reproducing the results

All scripts assume the working directory is the repository root (e.g. open `CointByParts.Rproj` in RStudio, which sets this automatically).

**Fastest path — just the empirical application**, using the calibration outputs already included in `output/cv_output/`:

```r
source("analysis/empirical_application.R")
```

**Full pipeline from scratch**, including re-running the calibration grid (slow — the full grid is 12 × 6 × 3 = 216 combinations × 10,000 replications):

```r
source("R/validate_fast_vs_r.R")      # sanity-check the C++ path first
source("R/run_generation_cv.R")       # full calibration grid (slow)
source("R/fit_response_surface.R")    # fit response surface from the grid
source("R/motivate_simulated_break.R")# motivating simulation figure
source("analysis/empirical_application.R")
```

## Dependencies

```r
install.packages(c(
  "data.table", "ggplot2", "urca", "Rcpp", "RcppArmadillo",
  "pbmcapply", "patchwork"
))
```

A working C++ compiler toolchain is required for `Rcpp::sourceCpp()`
(Rtools on Windows, Xcode Command Line Tools on macOS, `build-essential`
on Linux).

`.eps` figure output uses the `cairo_ps` graphics device, which requires
Cairo/X11 system libraries (on macOS, install XQuartz from
https://www.xquartz.org/ and restart afterward).

## Notes on intermediate output

`output/cv_output/checkpoints/` contains one file per `(T, n_regressors,
scan_type)` combination from the calibration grid, kept so the grid can be
inspected or resumed without recomputing from scratch. The raw
per-replication null-simulation draws (`critical_values_null_tmp/`) are
*not* included, as they are large, fully regenerable by `R/run_generation_cv.R`,
and not needed to reproduce any paper output.
