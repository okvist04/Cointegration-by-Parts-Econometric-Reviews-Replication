# =====================================================================
# run_size_power.R  --  empirical size and power evidence.
#
#   Rscript run_size_power.R
#
# SIZE. Rejection rates under the null at sample sizes NOT on the
# simulation grid, using critical values read off the fitted response
# surface -- the only check that the surface is usable at arbitrary T.
#
# POWER. Rejection rates against a break (both directions) and against
# an interior cointegrating window, benchmarked against the static
# full-sample Engle-Granger test and (optionally) Gregory-Hansen.
# Reported at nominal critical values from the response surface, and
# size-adjusted at critical values simulated at that exact T.
#
# Power is where Proposition 3.1 gets its evidence: FIEG should have no
# power against a reverse break, BIEG none against a forward break,
# GIEG power against both, and only GIEG any power against an interior
# window.
#
# Environment (defaults in brackets):
#   TABLES     [results/tables]  where make_tables.R wrote surface_main.csv
#   OUTDIR     [results/sizepower]
#   CASES      [c,ct]
#   T_SIZE     [90,140,350,600]      off-grid sample sizes
#   N_SIZE     [1,3]
#   T_POWER    [100,200,400]
#   NREP_SIZE  [5000]   NREP_POWER [2000]   NREP_NULL [5000]
#   R0 [0.15]  FLOOR_OBS [0]  PMAX_MULT [12]  DOF_RATIO [4]  BASE_SEED [20260907]
#   RHOS       [0.3,0.5,0.8]     TAUS [0.3,0.5,0.7]
#   GH         [1]     include the Gregory-Hansen benchmark (T <= 250)
#   MC_CORES   [parallel::detectCores() - 1]
# =====================================================================

source("R/engine.R")
source("R/response_surface.R")

ge_ <- function(n, d) Sys.getenv(n, d)
gi_ <- function(n, d) as.integer(ge_(n, d))
gf_ <- function(n, d) as.numeric(ge_(n, d))
ilist <- function(n, d) as.integer(strsplit(ge_(n, d), ",")[[1]])
flist <- function(n, d) as.numeric(strsplit(ge_(n, d), ",")[[1]])

TABLES     <- ge_("TABLES", "results/tables")
OUTDIR     <- ge_("OUTDIR", "results/sizepower")
CASES      <- trimws(strsplit(ge_("CASES", "c,ct"), ",")[[1]])
T_SIZE     <- ilist("T_SIZE", "90,140,350,600")
N_SIZE     <- ilist("N_SIZE", "1,3")
T_POWER    <- ilist("T_POWER", "100,200,400")
NREP_SIZE  <- gi_("NREP_SIZE", 5000)
NREP_POWER <- gi_("NREP_POWER", 2000)
NREP_NULL  <- gi_("NREP_NULL", 5000)
R0         <- gf_("R0", 0.15)
FLOOR_OBS  <- gi_("FLOOR_OBS", 0)
PMAX_MULT  <- gf_("PMAX_MULT", 12.0)
DOF_RATIO  <- gf_("DOF_RATIO", 4.0)
BASE_SEED  <- gi_("BASE_SEED", 20260907)
RHOS       <- flist("RHOS", "0.3,0.5,0.8")
TAUS       <- flist("TAUS", "0.3,0.5,0.7")
USE_GH     <- gi_("GH", 1) == 1
LEVELS     <- c(0.01, 0.05, 0.10)
MC_CORES   <- gi_("MC_CORES", max(1L, parallel::detectCores() - 1L))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

rejrate <- function(v, cv) {
  f <- v[is.finite(v)]
  if (length(f) == 0) return(NA_real_)
  mean(f < cv)
}

load_surfaces <- function() {
  f <- file.path(TABLES, "surface_main.csv")
  if (!file.exists(f)) stop(sprintf("%s not found -- run run_cv_grid.R then make_tables.R first", f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  S <- new.env()
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    key <- paste(r$case, r$N, r$test, round(r$prob, 4), sep = "|")
    S[[key]] <- list(beta = c(r$beta_inf, r$b1, r$b2, r$b3), se = c(r$se_inf, r$se1, r$se2, r$se3),
                      sigma = r$sigma, cubic = r$cubic, dof = r$dof, nobs = r$nobs)
  }
  S
}

nominal_cv <- function(S, case, N, test, level, T) {
  key <- paste(case, N, test, round(level, 4), sep = "|")
  if (is.null(S[[key]])) return(NA_real_)
  cv_at(S[[key]], T)
}

# One batch of replications under a given DGP kind. `kind` is "null",
# "break", or "window"; `params` carries kind-specific arguments.
run_batch <- function(kind, params, T, N, case, nrep, seed, with_gh = FALSE) {
  case_code <- case_to_code(case)
  m <- min_window(T, R0, FLOOR_OBS)
  set.seed(seed)
  fn <- if (with_gh) {
    switch(kind,
      null  = function(r) rep_null_gh_cpp(T, N, case_code, m, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO),
      break_ = function(r) rep_break_gh_cpp(T, N, case_code, m, params$tau0, params$direction, params$rho,
                                             pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO),
      window = function(r) rep_window_gh_cpp(T, N, case_code, m, params$tau1, params$tau2, params$rho,
                                              pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO)
    )
  } else {
    switch(kind,
      null  = function(r) rep_null_cpp(T, N, case_code, m, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO),
      break_ = function(r) rep_break_cpp(T, N, case_code, m, params$tau0, params$direction, params$rho,
                                          pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO),
      window = function(r) rep_window_cpp(T, N, case_code, m, params$tau1, params$tau2, params$rho,
                                           pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO)
    )
  }
  reps <- parallel_reps(fn, nrep, mc.cores = MC_CORES)
  fields <- c("fieg", "bieg", "gieg", "eg")
  if (with_gh) fields <- c(fields, "gh")
  collect_reps(reps, fields)
}

main <- function() {
  S <- load_surfaces()
  cat(strrep("=", 70), "\n")
  cat(sprintf("Size and power   cores %d   r0 %.2f floor %d\n", MC_CORES, R0, FLOOR_OBS))
  cat(strrep("=", 70), "\n")

  # ---------------- size --------------------------------------------
  size_rows <- list()
  for (case in CASES) for (T in T_SIZE) for (N in N_SIZE) {
    seed <- BASE_SEED + 31L * T + 17L * N + 5L
    t0 <- Sys.time()
    res <- run_batch("null", NULL, T, N, case, NREP_SIZE, seed)
    el <- as.numeric(Sys.time() - t0, units = "secs")
    for (test in c("FIEG", "BIEG", "GIEG")) for (lv in LEVELS) {
      cv <- nominal_cv(S, case, N, test, lv, T)
      size_rows[[length(size_rows) + 1]] <- data.frame(
        case = case, T = T, N = N, test = test, level = lv, cv_surface = cv,
        rejection = if (is.finite(cv)) rejrate(res[[tolower(test)]], cv) else NA_real_,
        nrep = NREP_SIZE
      )
    }
    saveRDS(list(stats = res, case = case, T = T, N = N, nrep = NREP_SIZE),
            file.path(OUTDIR, sprintf("size_case-%s_T-%d_N-%d.rds", case, T, N)))
    gieg5 <- rejrate(res$gieg, nominal_cv(S, case, N, "GIEG", 0.05, T))
    cat(sprintf("[size ] case %-2s T=%4d N=%d  %6.1f s  GIEG 5%% rejection %.3f\n", case, T, N, el, gieg5))
    flush.console()
    write.csv(do.call(rbind, size_rows), file.path(OUTDIR, "size.csv"), row.names = FALSE)
  }

  # ---------------- power ---------------------------------------------
  power_rows <- list()
  for (case in CASES) for (T in T_POWER) {
    N <- 1
    with_gh <- USE_GH && T <= 250
    t0 <- Sys.time()
    null_res <- run_batch("null", NULL, T, N, case, NREP_NULL, BASE_SEED + 991L * T, with_gh = with_gh)
    el <- as.numeric(Sys.time() - t0, units = "secs")
    matched <- list()
    for (t in names(null_res)) for (lv in LEVELS) {
      key <- paste(t, lv, sep = "|")
      matched[[key]] <- stats::quantile(null_res[[t]][is.finite(null_res[[t]])], lv, names = FALSE, type = 7)
    }
    saveRDS(list(stats = null_res, case = case, T = T, nrep = NREP_NULL),
            file.path(OUTDIR, sprintf("powernull_case-%s_T-%d.rds", case, T)))
    cat(sprintf("[null ] case %-2s T=%4d  %6.1f s  size-matched GIEG 5%% = %.3f\n",
                case, T, el, matched[["gieg|0.05"]])); flush.console()

    configs <- list()
    for (tau in TAUS) for (rho in RHOS) {
      configs[[length(configs) + 1]] <- list(name = "forward", tau = tau, rho = rho, direction = 0L)
      configs[[length(configs) + 1]] <- list(name = "reverse", tau = tau, rho = rho, direction = 1L)
    }
    for (rho in RHOS) {
      configs[[length(configs) + 1]] <- list(name = "window", tau = 0.5, rho = rho)
    }

    for (cfg in configs) {
      seed <- BASE_SEED + 7919L * T + round(1000 * cfg$tau) + round(100 * cfg$rho) +
        (if (cfg$name == "forward") 1L else if (cfg$name == "reverse") 2L else 3L) * 13L
      t1 <- Sys.time()
      if (cfg$name == "window") {
        res <- run_batch("window", list(tau1 = 0.25, tau2 = 0.75, rho = cfg$rho), T, N, case, NREP_POWER, seed,
                          with_gh = with_gh)
      } else {
        res <- run_batch("break_", list(tau0 = cfg$tau, direction = cfg$direction, rho = cfg$rho),
                          T, N, case, NREP_POWER, seed, with_gh = with_gh)
      }
      el2 <- as.numeric(Sys.time() - t1, units = "secs")
      for (test in names(res)) for (lv in LEVELS) {
        nom <- if (test %in% c("gh", "eg")) NA_real_ else nominal_cv(S, case, N, toupper(test), lv, T)
        matched_key <- paste(test, lv, sep = "|")
        matched_cv <- if (!is.null(matched[[matched_key]])) matched[[matched_key]] else NA_real_
        power_rows[[length(power_rows) + 1]] <- data.frame(
          case = case, T = T, N = N, dgp = cfg$name, tau = cfg$tau, rho = cfg$rho, test = toupper(test), level = lv,
          rej_nominal = if (is.finite(nom)) rejrate(res[[test]], nom) else NA_real_,
          rej_sizeadj = if (is.finite(matched_cv)) rejrate(res[[test]], matched_cv) else NA_real_,
          nrep = NREP_POWER
        )
      }
      saveRDS(list(stats = res, case = case, T = T, dgp = cfg$name, tau = cfg$tau, rho = cfg$rho),
              file.path(OUTDIR, sprintf("power_case-%s_T-%d_%s_tau-%.2f_rho-%.2f.rds",
                        case, T, cfg$name, cfg$tau, cfg$rho)))
      f5 <- function(t) rejrate(res[[t]], matched[[paste(t, "0.05", sep = "|")]])
      cat(sprintf("[power] case %-2s T=%4d %-7s tau=%.2f rho=%.2f  %6.1f s  size-adj 5%%: FIEG %.3f BIEG %.3f GIEG %.3f EG %.3f\n",
                  case, T, cfg$name, cfg$tau, cfg$rho, el2, f5("fieg"), f5("bieg"), f5("gieg"), f5("eg")))
      flush.console()
      write.csv(do.call(rbind, power_rows), file.path(OUTDIR, "power.csv"), row.names = FALSE)
    }
  }
  cat(sprintf("\nwrote size.csv and power.csv to %s\n", OUTDIR))
}

main()
