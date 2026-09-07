# ============================================================
# run_calibration.R
# ============================================================

library(data.table)
library(pbmcapply)

source("R/generate_cv.R")

T_GRID <- c(50, 75, 100, 150, 200, 250, 300, 400, 500, 700, 900, 1000)

N_REGRESSORS_GRID <- 1:6

N_REP <- 10000

SCAN_TYPES <- c("fw", "bw", "all")

PROBS <- c(0.01, 0.05, 0.10)

BASE_SEED <- 20260430

N_CORES <- max(1, parallel::detectCores() - 10)

OUTPUT_DIR <- "output/cv_output"

CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("====================================\n")
cat("Calibrating null critical values\n")
cat("====================================\n")
cat("T grid       :", T_GRID, "\n")
cat("N_regressors :", N_REGRESSORS_GRID, "\n")
cat("Replications :", N_REP, "\n")
cat("Cores        :", N_CORES, "\n")
cat("Scan types   :", SCAN_TYPES, "\n")
cat(
  "Combinations :",
  length(T_GRID) * length(N_REGRESSORS_GRID) * length(SCAN_TYPES),
  "\n"
)
cat("Output folder:", normalizePath(OUTPUT_DIR), "\n")
cat(
  "Note: each combination is checkpointed to '", CHECKPOINT_DIR,
  "/' independently. If this script is interrupted, just re-run it --\n",
  "completed combinations are loaded from disk instead of recomputed.\n",
  sep = ""
)

started_at <- Sys.time()

results <- calibrate_null_critical_values_parallel(
  T_grid = T_GRID,
  n_rep = N_REP,
  n_regressors_grid = N_REGRESSORS_GRID,
  scan_types = SCAN_TYPES,
  probs = PROBS,
  base_seed = BASE_SEED,
  n_cores = N_CORES,
  checkpoint_dir = CHECKPOINT_DIR,
  engine = "cpp"
)

fwrite(
  results,
  file.path(OUTPUT_DIR, "critical_values_null_full_grid.csv")
)

elapsed <- Sys.time() - started_at

cat("\n====================================\n")
cat("Finished\n")
cat("====================================\n")
cat("Rows written:", nrow(results), "\n")
cat("Elapsed     :", elapsed, "\n")
