# =====================================================================
# run_application.R  --  Section 5 recomputed, with p-values and the
# robustness checks the text promises.
#
#   Rscript run_application.R
#
# Produces, for both deterministic specifications:
#   * FIEG, BIEG, GIEG with their minimizing windows;
#   * response-surface critical values AND a p-value from the numerical
#     distribution function;
#   * robustness over r0 in {0.10, 0.15, 0.20}, BIC vs AIC, and with the
#     normalization swapped (the test is not invariant to which series
#     is on the left-hand side);
#   * ADF pre-tests on levels and differences of both series.
#
# NOTE: this script supersedes the simpler analysis/empirical_application.R
# from earlier in this project -- that script's FIEG/BIEG/GIEG/recursive
# break-search logic is still valid, but this version adds the full
# r0/BIC-AIC/normalization robustness grid and NDF-based p-values that
# script did not have.
#
# Environment:
#   SEALEVEL_CSV, TEMP_CSV     paths; each needs a year column and a value column
#   SEALEVEL_YEAR_COL [Time]   SEALEVEL_VAL_COL [GMSL (mm)]
#   TEMP_YEAR_COL [Time]       TEMP_VAL_COL [Anomaly (deg C)]
#   YEAR_MIN [1880]  YEAR_MAX [2019]
#   FLOOR_OBS [0]  PMAX_MULT [12]  DOF_RATIO [4]
#     -- these MUST match the settings the critical values were generated with
#   TABLES [results/tables]    OUTDIR [results/application]
# =====================================================================

source("R/engine.R")
source("R/response_surface.R")

ge_ <- function(n, d) Sys.getenv(n, d)
gi_ <- function(n, d) as.integer(ge_(n, d))

DATA_DIR      <- "data"  # matches the earlier analysis/ script's data/ folder
SEALEVEL_CSV  <- ge_("SEALEVEL_CSV", file.path(DATA_DIR, "CSIRO_Recons_gmsl_yr_2019.csv"))
TEMP_CSV      <- ge_("TEMP_CSV", file.path(DATA_DIR, "HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"))
SL_YEAR       <- ge_("SEALEVEL_YEAR_COL", "Time")
SL_VAL        <- ge_("SEALEVEL_VAL_COL", "GMSL (mm)")
TP_YEAR       <- ge_("TEMP_YEAR_COL", "Time")
TP_VAL        <- ge_("TEMP_VAL_COL", "Anomaly (deg C)")
FLOOR_OBS     <- gi_("FLOOR_OBS", 0)
PMAX_MULT     <- as.numeric(ge_("PMAX_MULT", "12"))
DOF_RATIO     <- as.numeric(ge_("DOF_RATIO", "4"))
YEAR_MIN      <- gi_("YEAR_MIN", 1880)
YEAR_MAX      <- gi_("YEAR_MAX", 2019)
TABLES        <- ge_("TABLES", "results/tables")
OUTDIR        <- ge_("OUTDIR", "results/application")
LEVELS        <- c(0.01, 0.05, 0.10)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_series <- function(path, yearcol, valcol, what) {
  if (!file.exists(path)) {
    stop(sprintf(paste0(
      "%s file not found: %s\n",
      "Set %s (and the column names) to point at the series used in Section 5."),
      what, path, if (what == "sea level") "SEALEVEL_CSV" else "TEMP_CSV"))
  }
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!(yearcol %in% names(df))) stop(sprintf("column '%s' not in %s; columns are %s", yearcol, path, paste(names(df), collapse = ", ")))
  if (!(valcol %in% names(df))) stop(sprintf("column '%s' not in %s; columns are %s", valcol, path, paste(names(df), collapse = ", ")))
  y <- as.integer(floor(as.numeric(df[[yearcol]])))
  v <- as.numeric(df[[valcol]])
  keep <- !is.na(v) & is.finite(v)
  data.frame(year = y[keep], value = v[keep])
}

surface_set <- function(df, case, N, test) {
  sub <- df[df$case == case & df$N == N & df$test == test, ]
  sub <- sub[order(sub$prob), ]
  probs <- sub$prob
  ss <- lapply(seq_len(nrow(sub)), function(i) {
    r <- sub[i, ]
    list(beta = c(r$beta_inf, r$b1, r$b2, r$b3), se = c(r$se_inf, r$se1, r$se2, r$se3),
         sigma = r$sigma, cubic = r$cubic, dof = r$dof, nobs = r$nobs)
  })
  list(probs = probs, surfaces = ss)
}

adf_series <- function(v, trend = FALSE, pmax_mult = PMAX_MULT) {
  adf_series_cpp(v, trend, pmax_mult)
}

main <- function() {
  sl <- read_series(SEALEVEL_CSV, SL_YEAR, SL_VAL, "sea level")
  tp <- read_series(TEMP_CSV, TP_YEAR, TP_VAL, "temperature")
  j <- merge(sl, tp, by = "year", suffixes = c("", "_1"))
  j <- j[j$year >= YEAR_MIN & j$year <= YEAR_MAX, ]
  j <- j[order(j$year), ]
  years <- j$year
  T <- nrow(j)
  cat(sprintf("sample: %d--%d, T = %d\n\n", years[1], years[length(years)], T))

  xv <- as.numeric(j$value)    # sea level  (regressor, per Schmith et al.)
  yv <- as.numeric(j$value_1)  # temperature (regressand)

  cat("integration-order pre-tests\n")
  for (pair in list(list("sea level", xv), list("temperature", yv))) {
    nm <- pair[[1]]; v <- pair[[2]]
    a1 <- adf_series(v, trend = FALSE); a2 <- adf_series(v, trend = TRUE)
    d1 <- adf_series(diff(v), trend = FALSE)
    cat(sprintf("  %-12s level ADF(c) %7.3f [p=%d]  ADF(ct) %7.3f [p=%d]   difference ADF(c) %7.3f [p=%d]\n",
                nm, a1$tstat, a1$lag, a2$tstat, a2$lag, d1$tstat, d1$lag))
  }
  cat("\n")

  sf <- file.path(TABLES, "surface_ndf.csv")
  surfaces_df <- if (file.exists(sf)) read.csv(sf, stringsAsFactors = FALSE) else NULL
  if (is.null(surfaces_df)) warning(sprintf("no surface_ndf.csv in %s -- statistics only, no critical values or p-values", TABLES))

  out_rows <- list()
  for (case in c("c", "ct")) for (r0 in c(0.10, 0.15, 0.20)) for (aic in c(FALSE, TRUE)) {
    for (norm in list(list("temp on sl", yv, xv), list("sl on temp", xv, yv))) {
      nrm <- norm[[1]]; yy <- norm[[2]]; xx <- norm[[3]]
      case_code <- case_to_code(case)
      m <- min_window(T, r0, FLOOR_OBS)
      res <- scan_all3_data_cpp(yy, matrix(xx, ncol = 1), case_code, m,
                                 pmax_mult = PMAX_MULT, dof_ratio = DOF_RATIO, aic = aic)
      for (test_info in list(list("FIEG", res$fieg, c(1, res$fe)),
                              list("BIEG", res$bieg, c(res$bs, T)),
                              list("GIEG", res$gieg, c(res$gs, res$ge)))) {
        test <- test_info[[1]]; stat <- test_info[[2]]; win <- test_info[[3]]
        cvs <- rep(NA_real_, 3); pv <- NA_real_
        if (!is.null(surfaces_df)) {
          ss <- surface_set(surfaces_df, case, 1, test)
          if (length(ss$surfaces) > 0) {
            for (i in seq_along(LEVELS)) {
              k <- which(abs(ss$probs - LEVELS[i]) < 1e-6)[1]
              cvs[i] <- if (is.na(k)) NA_real_ else cv_at(ss$surfaces[[k]], T)
            }
            pv <- pvalue_from_surface(ss$probs, ss$surfaces, T, stat)
          }
        }
        out_rows[[length(out_rows) + 1]] <- data.frame(
          case = case, r0 = r0, criterion = if (aic) "AIC" else "BIC", normalization = nrm,
          test = test, stat = stat, win_start = years[win[1]], win_end = years[win[2]],
          cv1 = cvs[1], cv5 = cvs[2], cv10 = cvs[3], pvalue = pv
        )
      }
    }
  }
  out <- do.call(rbind, out_rows)
  write.csv(out, file.path(OUTDIR, "application.csv"), row.names = FALSE)

  cat("headline specification (r0 = 0.15, BIC, temperature on sea level):\n")
  print(out[out$r0 == 0.15 & out$criterion == "BIC" & out$normalization == "temp on sl", ])

  # ---- remaining-segment check (Section 4.4) -------------------------
  # After locating the GIEG window, the natural next step is to remove
  # it and re-run the tests (in particular GIEG) on what remains, asking
  # whether a further episode of cointegration is present elsewhere.
  # This is only meaningful if the remaining segments themselves clear
  # the minimum sample size the calibration grid supports (T = 100,
  # per Section 4.3 -- shorter series over-reject). Report the segment
  # lengths and whether that floor is cleared, rather than silently
  # running (or silently skipping) the recursive test.
  MIN_SAMPLE_T <- as.integer(Sys.getenv("MIN_SAMPLE_T", "100"))
  headline <- out[out$r0 == 0.15 & out$criterion == "BIC" & out$normalization == "temp on sl" & out$test == "GIEG", ]
  if (nrow(headline) == 1) {
    win_start <- headline$win_start[1]; win_end <- headline$win_end[1]
    before_n <- win_start - years[1]
    after_n  <- years[length(years)] - win_end
    cat(sprintf("\nremaining-segment check (Section 4.4): GIEG window %d--%d\n", win_start, win_end))
    cat(sprintf("  segment before window: %d observations (%d--%d)\n",
                before_n, years[1], win_start - 1))
    cat(sprintf("  segment after window : %d observations (%d--%d)\n",
                after_n, win_end + 1, years[length(years)]))
    cat(sprintf("  minimum sample size supported by the calibration grid: T = %d\n", MIN_SAMPLE_T))
    if (before_n < MIN_SAMPLE_T && after_n < MIN_SAMPLE_T) {
      cat("  -> BOTH remaining segments fall below this floor; the recursive\n")
      cat("     break search is not run on this sample (see Section 4.4).\n")
    } else {
      cat("  -> at least one remaining segment clears the floor; a recursive\n")
      cat("     break search on that segment would be meaningful here.\n")
    }
  }

  cat(sprintf("\nwrote application.csv to %s\n", OUTDIR))
}

main()
