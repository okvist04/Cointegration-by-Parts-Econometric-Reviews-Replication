# Cointegration by Parts

Replication package for "Cointegration by Parts" (submitted to *Econometric
Reviews*, Special Issue in Honor of James G. MacKinnon).

This paper develops the **Cointegration by Parts** test — a recursive,
window-based procedure for detecting whether a cointegrating relationship
breaks down or is unstable over time, with finite-sample critical values
from a MacKinnon-style response-surface regression.

## The engine

`R/cbp_engine.cpp` is an Rcpp/C++ implementation of the recursive scan
(FIEG, BIEG, GIEG) across all three deterministic cases (no term /
constant / constant + trend). Every candidate ADF lag order is scored
on a common sample, and a degrees-of-freedom floor prevents short-window
numerical instability. `R/engine.R` wraps it with an R-level API
(DGPs, a parallel-replication helper, case-code mapping); `R/response_surface.R`
fits the response surfaces and numerical distribution function on top
of the calibration grid's output.

**Plain R would be far too slow for the actual calibration grid**
(50-200x slower than compiled code on tight nested loops like this),
which is why the hot path is C++, not an R script — everything else
(orchestration, tables, figures) is ordinary R.

**Validation**: the core engine was validated standalone (a portable
C++ test harness, no R/Rcpp needed) against a naive brute-force
reference before being wrapped for Rcpp — every window statistic across
all three deterministic cases and N = 1, 3 agreed to floating-point
precision (max relative error 1e-13 to 1e-15). `calibration/verify_engine.R`
reruns this same check from R and must print `ALL CHECKS PASSED` before
you trust anything downstream.

## Repository structure

```
R/
  cbp_engine.cpp, engine.R, response_surface.R   The engine (see above).
data/
  CSIRO_Recons_gmsl_yr_2019.csv
  HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv
calibration/
  verify_engine.R        RUN THIS FIRST -- correctness gate.
  run_cv_grid.R          The critical-value grids (long-running).
  make_tables.R          Quantiles, response surfaces, NDF, LaTeX.
  check_tail_fit.R       Does the grid reach its asymptote?
  check_grid_start.R     Is the smallest T distorting the surface?
  run_size_power.R       Empirical size (off-grid) and power.
  run_stage2_null.R      Valid critical values for the second-stage
                         subsample test.
  run_datestamp.R        How well is the cointegrating window located?
  run_application.R      Section 5, with the full r0 / BIC-AIC /
                         normalization robustness grid, NDF p-values,
                         and a remaining-segment check (Section 4.4)
                         reporting whether a recursive break search on
                         what's left of the sample would even be
                         meaningful, given the calibration grid's T = 100
                         floor.
paper/
  gen_tables.R, gen_surface.R, ndf_check.R      LaTeX table fragments.
  ndf_fig_data.R, ndf_fig.R                     Figure 2 (the NDF figure).
  series_fig.R                                  GMST_GMTA.eps: the raw
                                                GMSL/temperature series.
  windows_fig.R                                 cointegration_windows.eps:
                                                the application figure,
                                                window read directly from
                                                run_application.R's output.
results/                 Everything lands here, created on first run.
run_all.sh               The full pipeline, in order, resumable.
CointByParts.Rproj
```

## How the caching works

Every long-running `calibration/*.R` script writes one `.rds` file per
unit of work to `results/<stage>/` and skips it if that file already
exists. Once you've run the grid once, `make_tables.R`, every
`check_*.R` diagnostic, and everything in `paper/` only ever read
already-computed `.rds`/`.csv` files — changing a plot's styling or a
table's formatting never re-triggers a simulation. Only deleting the
relevant `results/` file (or folder) does.

## Running it

All scripts assume the working directory is the repository root.

**Correctness gate (run this first):**
```bash
Rscript calibration/verify_engine.R
```

**Fastest path to application results**, using whichever `results/`
already exists:
```bash
Rscript calibration/run_application.R
```

**Full pipeline from scratch:**
```bash
chmod +x run_all.sh
mkdir -p logs
nohup ./run_all.sh > logs/run_all.log 2>&1 &
tail -f logs/run_all.log
```

Override any setting via environment variables, e.g.:
```bash
R0=0.20 ./run_all.sh
CASES=c T_GRID=100,200,300 NREP=2000 Rscript calibration/run_cv_grid.R
```

**The two application figures** (`paper/series_fig.R` and
`paper/windows_fig.R`) both read `data/` directly with base R (the
`postscript` device — no extra packages needed). `windows_fig.R`
additionally reads `results/application/application.csv` for the
shaded window, so run `calibration/run_application.R` first.

## Dependencies

```r
install.packages(c("Rcpp", "parallel"))
```

That's the whole list — every script here uses only base R plus these
two packages. `Rcpp::sourceCpp()` (called by `R/engine.R`) needs a
working C++ compiler toolchain: Rtools on Windows, Xcode Command Line
Tools on macOS, `build-essential` on Linux. It compiles
`R/cbp_engine.cpp` once, on first `source("R/engine.R")` in a session.

Parallelism uses `parallel::mclapply()` (fork-based, Mac/Linux only —
on Windows every script automatically falls back to serial `lapply()`;
set `MC_CORES=1` anywhere to force serial execution on any platform).

## A couple of things worth knowing

- **Normal cdf/quantile**: `R/response_surface.R` uses R's built-in
  `pnorm()`/`qnorm()` rather than a hand-rolled rational approximation.
- **RNG / reproducibility**: uses R's own RNG via `set.seed()`.
  Reproducible given a fixed seed and a fixed `MC_CORES`, but not
  bitwise comparable across different `MC_CORES` values — a general
  limitation of fork-based parallel RNG, not specific to this code. For
  a fully deterministic single-machine record, run with `MC_CORES=1`.
