# =====================================================================
# run_cv_grid.R  --  null critical values for FIEG, BIEG and GIEG.
#
#   Rscript run_cv_grid.R
#
# Configured entirely through environment variables (defaults in brackets):
#
#   CASES      [c,ct]   deterministic cases of Equation (3.2) to run
#   T_GRID     [50,75,100,150,200,250,300,400,500,700,900,1000]
#   N_GRID     [1,2,3,4,5,6]
#   NREP       [10000]  Monte Carlo replications per cell
#   R0         [0.15]   minimum window fraction
#   FLOOR_OBS  [0]      hard floor on the window length; 0 = pure fraction
#   PMAX_MULT  [12]     Schwert constant in p_max = mult*(n/100)^(1/4)
#   BASE_SEED  [20260907]
#   OUTDIR     [results/cv]
#   DRIFT      [0]      drift in the null DGP (Case ct invariance check only)
#   DOF_RATIO  [0]      require effective sample >= ratio * (p+1) in the ADF
#                       regression; 0 keeps the paper's original rule
#   MC_CORES   [parallel::detectCores() - 1]
#
# Every cell writes ONE .rds file holding ALL of its simulated minima,
# not just three quantiles -- this is what makes the numerical
# distribution function, the Monte Carlo standard errors, and the
# bootstrap response-surface weights free afterwards in make_tables.R:
# they are read off the stored draws.
#
# Cells whose output file already exists are skipped, so the job is
# resumable: kill it and restart it and it picks up where it stopped.
# This is the "cache heavy steps to disk" behavior the whole R port was
# built around -- once a cell's .rds exists, make_tables.R and every
# downstream script never re-run the simulation to get a plot.
# =====================================================================

source("R/engine.R")

getenv <- function(name, default) Sys.getenv(name, default)
getint <- function(name, default) as.integer(getenv(name, default))
getflt <- function(name, default) as.numeric(getenv(name, default))
getlist_int <- function(name, default) as.integer(strsplit(getenv(name, default), ",")[[1]])

CASES     <- trimws(strsplit(getenv("CASES", "c,ct"), ",")[[1]])
T_GRID    <- getlist_int("T_GRID", "50,75,100,150,200,250,300,400,500,700,900,1000")
N_GRID    <- getlist_int("N_GRID", "1,2,3,4,5,6")
NREP      <- getint("NREP", 10000)
R0        <- getflt("R0", 0.15)
FLOOR_OBS <- getint("FLOOR_OBS", 0)
PMAX_MULT <- getflt("PMAX_MULT", 12.0)
BASE_SEED <- getint("BASE_SEED", 20260907)
OUTDIR    <- getenv("OUTDIR", "results/cv")
DRIFT     <- getflt("DRIFT", 0.0)
DOF_RATIO <- getflt("DOF_RATIO", 4.0)
MC_CORES  <- getint("MC_CORES", max(1L, parallel::detectCores() - 1L))
PROBS     <- c(0.01, 0.05, 0.10)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cellfile <- function(case, T, N) {
  suffix <- paste0(
    if (DRIFT != 0) "_drift" else "",
    if (DOF_RATIO != 0) sprintf("_dof-%.1f", DOF_RATIO) else ""
  )
  file.path(OUTDIR, sprintf("minima_case-%s_T-%d_N-%d_r0-%.2f_floor-%d_R-%d%s.rds",
                             case, T, N, R0, FLOOR_OBS, NREP, suffix))
}

# One cell of the grid: NREP null replications at (case, T, N).
run_cell <- function(case, T, N) {
  m <- min_window(T, R0, FLOOR_OBS)
  case_code <- case_to_code(case)
  cellseed <- BASE_SEED + 1000003L * T + 7919L * N +
    104729L * (if (case == "n") 1L else if (case == "c") 2L else 3L)
  set.seed(cellseed)
  reps <- parallel_reps(function(r) {
    rep_null_cpp(T, N, case_code, m, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO, drift = DRIFT)
  }, NREP, mc.cores = MC_CORES)
  co <- collect_reps(reps, c("fieg", "bieg", "gieg", "fe", "bs", "gs", "ge"))
  list(fieg = co$fieg, bieg = co$bieg, gieg = co$gieg,
       fe = co$fe, bs = co$bs, gs = co$gs, ge = co$ge, min_obs = m)
}

main <- function() {
  cat(strrep("=", 70), "\n")
  cat("Null critical values for FIEG / BIEG / GIEG\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("  cases      : %s\n", paste(CASES, collapse = ", ")))
  cat(sprintf("  T grid     : %s\n", paste(T_GRID, collapse = ", ")))
  cat(sprintf("  N grid     : %s\n", paste(N_GRID, collapse = ", ")))
  cat(sprintf("  reps       : %d\n", NREP))
  cat(sprintf("  r0 / floor : %.2f / %d observations\n", R0, FLOOR_OBS))
  cat(sprintf("  p_max rule : floor(%.1f*(n/100)^0.25), n = window length\n", PMAX_MULT))
  cat(sprintf("  cores      : %d\n", MC_CORES))
  cat(sprintf("  output     : %s\n", OUTDIR))
  if (DOF_RATIO != 0) cat(sprintf("  DOF_RATIO  : %.1f\n", DOF_RATIO))
  if (DRIFT != 0) cat(sprintf("  DRIFT      : %.3f  (invariance check run)\n", DRIFT))
  cat(strrep("=", 70), "\n")

  summary_rows <- list()
  t0 <- Sys.time()
  for (case in CASES) for (T in T_GRID) for (N in N_GRID) {
    f <- cellfile(case, T, N)
    if (file.exists(f)) {
      cat(sprintf("[skip] case %-2s T=%4d N=%d  (exists)\n", case, T, N)); flush.console()
      next
    }
    t1 <- Sys.time()
    cell <- run_cell(case, T, N)
    el <- as.numeric(Sys.time() - t1, units = "secs")
    saveRDS(list(fieg = cell$fieg, bieg = cell$bieg, gieg = cell$gieg,
                 fe = cell$fe, bs = cell$bs, gs = cell$gs, ge = cell$ge,
                 case = case, T = T, N = N, nrep = NREP, r0 = R0,
                 floor_obs = FLOOR_OBS, min_obs = cell$min_obs,
                 pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO,
                 base_seed = BASE_SEED, drift = DRIFT,
                 common_sample_bic = TRUE, seconds = el, finished = as.character(Sys.time())),
            f)
    for (nm in c("FIEG", "BIEG", "GIEG")) {
      v <- switch(nm, FIEG = cell$fieg, BIEG = cell$bieg, GIEG = cell$gieg)
      fin <- v[is.finite(v)]
      for (p in PROBS) {
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          case = case, T = T, N = N, test = nm, prob = p,
          cv = stats::quantile(fin, p, names = FALSE, type = 7),
          n_finite = length(fin), min_obs = cell$min_obs, seconds = el
        )
      }
    }
    gieg5 <- stats::quantile(cell$gieg[is.finite(cell$gieg)], 0.05, names = FALSE, type = 7)
    cat(sprintf("[done] case %-2s T=%4d N=%d  %8.1f s  GIEG 5%% = %7.3f\n", case, T, N, el, gieg5))
    flush.console()
    summary_df <- do.call(rbind, summary_rows)
    write.csv(summary_df, file.path(OUTDIR, sprintf("summary_case-%s_N-%s_R-%d.csv",
              paste(CASES, collapse = "-"), paste(N_GRID, collapse = "-"), NREP)), row.names = FALSE)
  }
  cat(sprintf("\nTotal wall time: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  cat("Next: Rscript make_tables.R\n")
}

main()
