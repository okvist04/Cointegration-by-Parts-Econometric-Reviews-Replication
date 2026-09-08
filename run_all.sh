#!/usr/bin/env bash
# =====================================================================
# run_all.sh -- everything the calibration pipeline needs, in order.
#
#   cd "/path/to/this/repo"
#   chmod +x run_all.sh
#   nohup ./run_all.sh > logs/run_all.log 2>&1 &
#   tail -f logs/run_all.log
#
# Safe to kill and restart: run_cv_grid.R skips cells that already have
# an output .rds file, and every step writes into results/ under its
# own name. This resumability is the whole point of the R port's
# caching design: once a step's output exists, no downstream script
# (tables, figures, robustness checks) ever re-runs the simulation to
# regenerate a plot.
# =====================================================================
set -u
cd "$(dirname "$0")"
mkdir -p logs results

# ---- shared settings: every step below uses the same statistic --------
export R0=${R0:-0.15}                 # minimum window FRACTION (never a fixed count)
export FLOOR_OBS=${FLOOR_OBS:-0}      # no floor: a floor makes r0 vary with T
export PMAX_MULT=${PMAX_MULT:-12}     # Schwert constant in the ADF lag bound
export DOF_RATIO=${DOF_RATIO:-4}      # ADF keeps neff >= 4(p+1)
export T_GRID=${T_GRID:-100,150,200,250,300,400,500,700,900,1000}
export BASE_SEED=${BASE_SEED:-20260907}
export NREP=${NREP:-10000}
MC_CORES=${MC_CORES:-$(Rscript -e 'cat(max(1L, parallel::detectCores() - 1L))' 2>/dev/null || echo 4)}
export MC_CORES

echo "=== settings: r0=$R0 floor=$FLOOR_OBS pmax_mult=$PMAX_MULT dof_ratio=$DOF_RATIO"
echo "=== T grid: $T_GRID"
echo "=== cores=$MC_CORES reps=$NREP  started $(date)"

step () {                              # step <name> <command...>
  local name=$1; shift
  local log="logs/${name}.log"
  if [ -f "logs/${name}.done" ]; then echo "--- skip $name (done)"; return; fi
  echo "--- $name  ($(date +%H:%M))"
  if "$@" > "$log" 2>&1; then
    touch "logs/${name}.done"; echo "    ok -> $log"
  else
    echo "    FAILED -- see $log"; tail -20 "$log"
    [ "${SOFT:-0}" = "1" ] && { echo "    (optional step, continuing)"; return 0; }
    exit 1
  fi
}

# STEP 0: the engine must agree with a naive brute-force reference.
step 00_verify Rscript calibration/verify_engine.R

# STEP 1: the critical value grids. Case c and Case ct. Long-running.
step 01_cv_case_c  env CASES=c  Rscript calibration/run_cv_grid.R
step 02_cv_case_ct env CASES=ct Rscript calibration/run_cv_grid.R

# STEP 2: quantiles, Monte Carlo standard errors, response surfaces
# with standard errors, and the numerical distribution function.
step 03_tables Rscript calibration/make_tables.R

# STEP 3: empirical size off the grid, and power against both break
# directions and an interior window.
step 04_size_power Rscript calibration/run_size_power.R

# STEP 4: valid critical values for the second-stage subsample test, at
# the two sample sizes the application uses.
step 05_stage2_T140 env T=140 Rscript calibration/run_stage2_null.R
step 06_stage2_T96  env T=96  Rscript calibration/run_stage2_null.R

# STEP 5: how well the minimizing window locates the true one.
step 07_datestamp Rscript calibration/run_datestamp.R

# STEP 6: does the GIEG column flatten? N = 1 only, fewer replications.
step 08_gieg_extension env CASES=c,ct T_GRID=1500,2000 N_GRID=1 NREP=2500 \
  Rscript calibration/run_cv_grid.R

# STEP 7: robustness diagnostics on the fitted surface (optional; only
# meaningful once the full grid above has finished).
SOFT=1 step 09_tail_fit    Rscript calibration/check_tail_fit.R
SOFT=1 step 10_grid_start  Rscript calibration/check_grid_start.R

# STEP 8: the application. Needs data/ (already included in this repo).
SOFT=1 step 11_application Rscript calibration/run_application.R

# STEP 9: paper tables and figures, from results/ written above.
SOFT=1 step 12_paper_tables    Rscript paper/gen_tables.R
SOFT=1 step 13_paper_surface   Rscript paper/gen_surface.R
SOFT=1 step 14_paper_ndf_check Rscript paper/ndf_check.R
SOFT=1 step 15_paper_ndf_data  Rscript paper/ndf_fig_data.R
SOFT=1 step 16_paper_ndf_fig   Rscript paper/ndf_fig.R
SOFT=1 step 17_paper_series    Rscript paper/series_fig.R
SOFT=1 step 18_paper_windows   Rscript paper/windows_fig.R

echo "=== finished $(date)"
echo "results/ now holds everything; paper/ holds the .tex fragments and figures."
