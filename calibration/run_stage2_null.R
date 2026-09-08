# =====================================================================
# run_stage2_null.R  --  valid critical values for the second-stage
# subsample test of Section 5.
#
#   Rscript run_stage2_null.R
#
# Table 5.2 re-runs the three tests on a subsample chosen BECAUSE the
# first-stage GIEG window excluded it, then compares the result to
# critical values simulated for an UNCONDITIONAL null. Selecting the
# sample on the outcome of the first test and then using unconditional
# critical values does not have the nominal level.
#
# This script simulates the conditional null directly: repeat the whole
# two-stage procedure on null data and take the quantiles of the
# second-stage statistic. Those are the critical values Table 5.2
# should be compared against.
#
# Environment:
#   T [140]  CASE [c]  NREP [10000]  R0 [0.15]  FLOOR_OBS [0]  DOF_RATIO [4]
#   MIN_SUB [30]   shortest subsample that is still tested
#   STAGE1_LEVEL [0.01]   level at which the first stage must reject
#   TABLES [results/tables]   OUTDIR [results/stage2]
#   MC_CORES [parallel::detectCores() - 1]
# =====================================================================

source("R/engine.R")
source("R/response_surface.R")

ge_ <- function(n, d) Sys.getenv(n, d)
gi_ <- function(n, d) as.integer(ge_(n, d))
gf_ <- function(n, d) as.numeric(ge_(n, d))

T         <- gi_("T", 140)
CASE      <- ge_("CASE", "c")
NREP      <- gi_("NREP", 10000)
R0        <- gf_("R0", 0.15)
FLOOR_OBS <- gi_("FLOOR_OBS", 0)
MIN_SUB   <- gi_("MIN_SUB", 30)
S1LEVEL   <- gf_("STAGE1_LEVEL", 0.01)
PMAX_MULT <- gf_("PMAX_MULT", 12.0)
DOF_RATIO <- gf_("DOF_RATIO", 4.0)
SEED      <- gi_("BASE_SEED", 20260907)
TABLES    <- ge_("TABLES", "results/tables")
OUTDIR    <- ge_("OUTDIR", "results/stage2")
LEVELS    <- c(0.01, 0.05, 0.10)
MC_CORES  <- gi_("MC_CORES", max(1L, parallel::detectCores() - 1L))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

surfaces <- function() {
  f <- file.path(TABLES, "surface_main.csv")
  if (!file.exists(f)) stop(sprintf("%s not found -- run make_tables.R first", f))
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

main <- function() {
  S <- surfaces()
  cv1 <- cv_at(S[[paste(CASE, 1, "GIEG", round(S1LEVEL, 4), sep = "|")]], T)
  cat(sprintf("two-stage conditional null: T = %d, case %s, stage-1 GIEG %.0f%% cv = %.3f\n",
              T, CASE, 100 * S1LEVEL, cv1))

  case_code <- case_to_code(CASE)
  m <- min_window(T, R0, FLOOR_OBS)

  # Design note: computing the first-stage scan and the second-stage scan
  # on a subsample both need the SAME underlying (y, x) draw, so each
  # replication draws once via dgp_null_cpp() and reuses it for both scans
  # (unlike rep_null_cpp(), which is a fast path that never returns y/x).
  set.seed(SEED + 13L * T)
  do_one <- function(r) {
    d <- dgp_null_cpp(T, 1)
    res <- scan_all3_data_cpp(d$y, d$x, case_code, m, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO)
    stat1 <- res$gieg
    if (!is.finite(stat1) || stat1 >= cv1) {
      return(list(stat1 = stat1, rej = FALSE, n2 = 0L, side = 0L,
                  s2_fieg = NA_real_, s2_bieg = NA_real_, s2_gieg = NA_real_))
    }
    gs <- res$gs; ge <- res$ge
    cands <- list(c(1, gs - 1, 1L), c(ge + 1, T, 2L))
    best <- 0L; bl <- 0L; bside <- 0L
    for (cc in cands) {
      n <- cc[2] - cc[1] + 1
      if (n >= MIN_SUB && n > bl) { best <- cc[1]; bl <- n; bside <- cc[3] }
    }
    if (bl == 0L) {
      return(list(stat1 = stat1, rej = TRUE, n2 = 0L, side = 0L,
                  s2_fieg = NA_real_, s2_bieg = NA_real_, s2_gieg = NA_real_))
    }
    idx <- best:(best + bl - 1)
    y2 <- d$y[idx]; x2 <- d$x[idx, , drop = FALSE]
    m2 <- min_window(bl, R0, FLOOR_OBS)
    if (m2 >= bl) {
      return(list(stat1 = stat1, rej = TRUE, n2 = bl, side = bside,
                  s2_fieg = NA_real_, s2_bieg = NA_real_, s2_gieg = NA_real_))
    }
    r2 <- scan_all3_data_cpp(y2, x2, case_code, m2, pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO)
    list(stat1 = stat1, rej = TRUE, n2 = bl, side = bside,
         s2_fieg = r2$fieg, s2_bieg = r2$bieg, s2_gieg = r2$gieg)
  }

  t0 <- Sys.time()
  reps <- parallel_reps(do_one, NREP, mc.cores = MC_CORES)
  el <- as.numeric(Sys.time() - t0, units = "secs")

  co <- collect_reps(reps, c("stat1", "n2", "side", "s2_fieg", "s2_bieg", "s2_gieg"))
  rej <- vapply(reps, function(r) isTRUE(r$rej), logical(1))
  keep <- (co$n2 > 0) & is.finite(co$s2_gieg)

  cat(sprintf("  %.1f s;  stage-1 rejects in %.3f of draws;  a testable subsample in %.3f\n",
              el, mean(rej), mean(keep)))
  cat(sprintf("  subsample length: median %s, mean %.1f\n",
              if (any(keep)) stats::median(co$n2[keep]) else 0, if (any(keep)) mean(co$n2[keep]) else NaN))

  nmed <- if (any(keep)) round(stats::median(co$n2[keep])) else 0
  out_rows <- list()
  for (nm in c("FIEG", "BIEG", "GIEG")) {
    v <- co[[paste0("s2_", tolower(nm))]][keep]
    v <- v[is.finite(v)]
    for (lv in LEVELS) {
      key <- paste(CASE, 1, nm, round(lv, 4), sep = "|")
      uncond <- if (nmed > 0 && !is.null(S[[key]])) cv_at(S[[key]], nmed) else NA_real_
      out_rows[[length(out_rows) + 1]] <- data.frame(
        test = nm, level = lv,
        conditional_cv = if (length(v) > 0) stats::quantile(v, lv, names = FALSE, type = 7) else NA_real_,
        unconditional_cv_at_median_n = uncond, n_used = length(v)
      )
    }
  }
  df <- do.call(rbind, out_rows)
  write.csv(df, file.path(OUTDIR, sprintf("stage2_case-%s_T-%d.csv", CASE, T)), row.names = FALSE)
  saveRDS(list(stat1 = co$stat1, rej = rej, n2 = co$n2, side = co$side,
               s2_fieg = co$s2_fieg, s2_bieg = co$s2_bieg, s2_gieg = co$s2_gieg,
               T = T, case = CASE, nrep = NREP, cv1 = cv1,
               r0 = R0, floor_obs = FLOOR_OBS, min_sub = MIN_SUB, stage1_level = S1LEVEL),
          file.path(OUTDIR, sprintf("stage2_case-%s_T-%d.rds", CASE, T)))
  cat("\nconditional vs unconditional critical values (second stage):\n")
  print(df)
  cat(sprintf("\nwrote %s\n", OUTDIR))
}

main()
