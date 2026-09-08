# =====================================================================
# run_datestamp.R  --  how well is the cointegrating window located?
#
#   Rscript run_datestamp.R
#
# Simulates data whose cointegrating window is known, then records
# where the GIEG minimum actually falls. Reported per configuration:
#   bias and RMSE of both endpoints, in sample fractions
#   P(minimizing window inside the true window)
#   overlap (intersection over union) with the true window
#
# Environment:
#   T_GRID [100,140,200,400]  CASES [c,ct]  RHOS [0.3,0.5,0.8]
#   TAU1 [0.25]  TAU2 [0.75]  NREP [2000]
#   R0 [0.15]  FLOOR_OBS [0]  OUTDIR [results/datestamp]
#   MC_CORES [parallel::detectCores() - 1]
# =====================================================================

source("R/engine.R")

ge_ <- function(n, d) Sys.getenv(n, d)
gi_ <- function(n, d) as.integer(ge_(n, d))
gf_ <- function(n, d) as.numeric(ge_(n, d))
ilist <- function(n, d) as.integer(strsplit(ge_(n, d), ",")[[1]])
flist <- function(n, d) as.numeric(strsplit(ge_(n, d), ",")[[1]])

T_GRID    <- ilist("T_GRID", "100,140,200,400")
CASES     <- trimws(strsplit(ge_("CASES", "c,ct"), ",")[[1]])
RHOS      <- flist("RHOS", "0.3,0.5,0.8")
TAU1      <- gf_("TAU1", 0.25)
TAU2      <- gf_("TAU2", 0.75)
NREP      <- gi_("NREP", 2000)
R0        <- gf_("R0", 0.15)
FLOOR_OBS <- gi_("FLOOR_OBS", 0)
PMAX_MULT <- gf_("PMAX_MULT", 12.0)
DOF_RATIO <- gf_("DOF_RATIO", 4.0)
SEED      <- gi_("BASE_SEED", 20260907)
OUTDIR    <- ge_("OUTDIR", "results/datestamp")
MC_CORES  <- gi_("MC_CORES", max(1L, parallel::detectCores() - 1L))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

main <- function() {
  rows <- list()
  for (case in CASES) for (T in T_GRID) for (rho in RHOS) {
    case_code <- case_to_code(case)
    m <- min_window(T, R0, FLOOR_OBS)
    seed <- SEED + 977L * T + round(100 * rho)
    set.seed(seed)
    t0 <- Sys.time()
    reps <- parallel_reps(function(r) {
      rep_window_cpp(T, 1, case_code, m, TAU1, TAU2, rho, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO)
    }, NREP, mc.cores = MC_CORES)
    el <- as.numeric(Sys.time() - t0, units = "secs")

    co <- collect_reps(reps, c("gs", "ge", "gieg", "t1", "t2"))
    t1 <- co$t1[1]; t2 <- co$t2[1]  # same for every rep (deterministic from T, TAU1, TAU2)
    ok <- (co$gs > 0) & (co$ge > 0)
    s <- co$gs[ok] / T; e <- co$ge[ok] / T
    inside <- mean(co$gs[ok] >= t1 + 1 & co$ge[ok] <= t2)
    idx_ok <- which(ok)
    ov <- vapply(idx_ok, function(i) {
      inter <- max(0, min(co$ge[i], t2) - max(co$gs[i], t1 + 1) + 1)
      uni <- max(co$ge[i], t2) - min(co$gs[i], t1 + 1) + 1
      inter / uni
    }, numeric(1))

    rows[[length(rows) + 1]] <- data.frame(
      case = case, T = T, rho = rho, tau1 = TAU1, tau2 = TAU2,
      bias_start = mean(s) - TAU1, rmse_start = sqrt(mean((s - TAU1)^2)),
      bias_end = mean(e) - TAU2, rmse_end = sqrt(mean((e - TAU2)^2)),
      p_inside = inside, mean_overlap = mean(ov),
      mean_len = mean(co$ge[ok] - co$gs[ok] + 1), nrep = NREP
    )
    saveRDS(list(gs = co$gs, ge = co$ge, gieg = co$gieg, T = T, case = case, rho = rho,
                 tau1 = TAU1, tau2 = TAU2, nrep = NREP, min_obs = m),
            file.path(OUTDIR, sprintf("datestamp_case-%s_T-%d_rho-%.2f.rds", case, T, rho)))
    cat(sprintf("[stamp] case %-2s T=%4d rho=%.2f  %6.1f s  inside %.3f  overlap %.3f  RMSE (%.3f, %.3f)\n",
                case, T, rho, el, inside, mean(ov), sqrt(mean((s - TAU1)^2)), sqrt(mean((e - TAU2)^2))))
    flush.console()
    write.csv(do.call(rbind, rows), file.path(OUTDIR, "datestamp.csv"), row.names = FALSE)
  }
  cat(sprintf("\nwrote datestamp.csv to %s\n", OUTDIR))
}

main()
