# ============================================================
# prepare_and_test_sealevel_temperature.R
#
# 1. Loads CSIRO global mean sea level (GMSL) data
# 2. Loads HadCRUT5 global temperature data and rebases it from its
#    native 1961-1990 reference period to 1850-1900
# 3. Aligns both series to a common annual sample
# 4. Runs FIEG(r0)/BIEG(r0,1)/GIEG(r0) on the result via run_scan_test.R
#
# I haven't seen your actual downloaded files, so the column-name
# detection below is defensive rather than assumed -- it prints what
# it finds and stops with a clear message if it can't match something,
# rather than silently guessing wrong.
# ============================================================

library(data.table); library(ggplot2); library(urca)

# ------------------------------------------------------------
# CONFIGURATION -- edit these to match your actual downloaded files
# ------------------------------------------------------------

CSIRO_FILE   <- "data/CSIRO_Recons_gmsl_yr_2019.csv"                              # <- set to your actual file
HADCRUT_FILE <- "data/HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"  # <- set to your actual file

# Which series is the dependent variable (y) and which is the
# regressor (X)? Sea level as y, temperature as X is the common
# convention in the semi-empirical sea-level literature (e.g.
# Rahmstorf 2007). Flip this if you'd rather match a specific
# specification (e.g. Schmith et al.'s reported causal direction).
Y_SERIES <- "temperature"   # "sea_level" or "temperature"

# Restrict the sample further if you want a specific window -- leave
# as -Inf/Inf to use the full overlap between the two files (printed
# below once both are loaded).
YEAR_MIN <- -Inf
YEAR_MAX <- Inf

R0 <- 0.15            # trimming fraction -- MUST match whatever your

# ------------------------------------------------------------
# 1. Load CSIRO global mean sea level
# ------------------------------------------------------------

gmsl_raw <- fread(CSIRO_FILE)
cat("CSIRO file columns:", paste(names(gmsl_raw), collapse = ", "), "\n")

time_col_gmsl <- grep("time|year|date", names(gmsl_raw), ignore.case = TRUE, value = TRUE)[1]
gmsl_col      <- grep("gmsl|sea.?level|adjusted", names(gmsl_raw), ignore.case = TRUE, value = TRUE)[1]

if (is.na(time_col_gmsl) || is.na(gmsl_col)) {
  stop(
    "Could not auto-detect the time/GMSL columns in the CSIRO file. ",
    "Set time_col_gmsl and gmsl_col manually to match the column names printed above."
  )
}
cat(sprintf("Using '%s' as time and '%s' as GMSL.\n", time_col_gmsl, gmsl_col))

gmsl <- gmsl_raw[, .(
  year = {
    tc <- get(time_col_gmsl)
    if (is.numeric(tc)) floor(tc) else as.integer(substr(as.character(tc), 1, 4))
  },
  sea_level = as.numeric(get(gmsl_col))
)]

# Collapse to one value per year (averages sub-annual timestamps, if any).
gmsl <- gmsl[, .(sea_level = mean(sea_level, na.rm = TRUE)), by = year]

cat(sprintf(
  "Sea level: %d annual observations, %d-%d\n",
  nrow(gmsl), min(gmsl$year), max(gmsl$year)
))

# ------------------------------------------------------------
# 2. Load HadCRUT5 and rebase to 1850-1900
# ------------------------------------------------------------

hadcrut_raw <- fread(HADCRUT_FILE)
cat("HadCRUT file columns:", paste(names(hadcrut_raw), collapse = ", "), "\n")

setnames(
  hadcrut_raw,
  old = c("Time", "Anomaly (deg C)"),
  new = c("time", "anomaly_1961_1990"),
  skip_absent = TRUE
)

if (!all(c("time", "anomaly_1961_1990") %in% names(hadcrut_raw))) {
  stop(
    "Could not find the expected 'Time' / 'Anomaly (deg C)' columns in the HadCRUT file. ",
    "Check its actual header names (printed above) and adjust the setnames() call."
  )
}

hadcrut_raw[, year := if (is.numeric(time)) as.integer(time) else as.integer(substr(as.character(time), 1, 4))]

baseline <- hadcrut_raw[year >= 1850 & year <= 1900]
if (nrow(baseline) == 0) {
  stop("No rows found for 1850-1900 in the HadCRUT file -- check that 'year' parsed correctly.")
}

offset <- mean(baseline$anomaly_1961_1990, na.rm = TRUE)
cat(sprintf(
  "HadCRUT5 offset applied (1961-1990 -> 1850-1900 baseline): %.4f degC\n", offset
))

hadcrut_raw[, anomaly_1850_1900 := anomaly_1961_1990 - offset]

# Collapse to one value per year (averages monthly values, if the file
# was monthly rather than annual).
temperature <- hadcrut_raw[, .(temperature = mean(anomaly_1850_1900, na.rm = TRUE)), by = year]

cat(sprintf(
  "Temperature: %d annual observations, %d-%d\n",
  nrow(temperature), min(temperature$year), max(temperature$year)
))

# ------------------------------------------------------------
# 3. Align on common years
# ------------------------------------------------------------

combined <- merge(gmsl, temperature, by = "year")
combined <- combined[year >= YEAR_MIN & year <= YEAR_MAX]
combined <- combined[is.finite(sea_level) & is.finite(temperature)]
setorder(combined, year)

if (nrow(combined) < 20) {
  stop(sprintf(
    "Only %d overlapping years after alignment -- check the two files actually cover comparable periods.",
    nrow(combined)
  ))
}

cat(sprintf(
  "\nCombined sample: %d annual observations, %d-%d\n",
  nrow(combined), min(combined$year), max(combined$year)
))

# ------------------------------------------------------------
# 4. Plot the aligned series
# ------------------------------------------------------------
#
# Sea level (mm) and temperature anomaly (degC) are on different
# scales, so to show them overlapping in one panel, temperature is
# linearly rescaled onto sea level's range purely for plotting -- the
# right-hand axis then back-transforms to show temperature in its own
# real units. This is the standard dual-axis idiom, but worth being
# upfront about in the paper if this figure is used there: the
# vertical alignment is chosen so the two RANGES match, not because of
# any principled unit conversion between them, so the apparent
# closeness of the two lines shouldn't be read as more meaningful than
# that -- it's a visualization choice, not a result.

range_sl   <- range(combined$sea_level, na.rm = TRUE)
range_temp <- range(combined$temperature, na.rm = TRUE)

rescale_lin <- function(x, from, to) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

combined[, temperature_scaled := rescale_lin(temperature, range_temp, range_sl)]

p_overlay <- ggplot(combined, aes(x = year)) +
  geom_line(aes(y = sea_level, colour = "Sea level (CSIRO GMSL)")) +
  geom_line(aes(y = temperature_scaled, colour = "Temperature anomaly (HadCRUT5)")) +
  scale_y_continuous(
    name = "Sea level (mm)",
    sec.axis = sec_axis(
      transform = ~ rescale_lin(., range_sl, range_temp),
      name = "Temperature anomaly (degC, 1850-1900 baseline)"
    )
  ) +
  scale_colour_manual(name = NULL, values = c("steelblue", "firebrick")) +
  labs(
    x = "Year",
    title = "Global mean sea level and land-ocean temperature anomaly",
    subtitle = sprintf("%d-%d, temperature rescaled onto sea level's range for display", min(combined$year), max(combined$year))
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_overlay)

ggsave("output/figures/GMST_GMTA.eps",
       p_overlay,
       device = cairo_ps,
       width = 8,
       height = 6,
       units = "in")

# ------------------------------------------------------------
# 4A. Run ADF tests on the full sample
# ------------------------------------------------------------
# Schmith, Johansen & Thejll (2012, J. Climate) report a univariate
# ADF test with a constant term only (no linear trend) and two FIXED
# lags of the differenced series (not selected via an information
# criterion), obtaining p = 0.83 for temperature and p = 0.99 for sea
# level -- in both cases, a clear failure to reject the unit-root
# null. This corresponds to urca::ur.df()'s type = "drift" model with
# lags = 2 fixed, we however do not restrict the lags to 2, we let BIC 
# choose.
#
# Note: urca::ur.df() does not report a p-value directly, only the
# test statistic against Dickey-Fuller critical values (as above for
# the drift/BIC case) -- so this replicates Schmith et al.'s exact
# specification and substantive conclusion (fail to reject the unit
# root for both series), but the specific p-values they report likely
# come from a different software implementation's interpolation and
# are not reproduced exactly here.

adf_gmsl <- ur.df(combined$sea_level, type = "drift", selectlags = "BIC")
adf_temp <- ur.df(combined$temperature, type = "drift", selectlags = "BIC")

summary(adf_gmsl) # test-statistics are 1.9111, 13.1937 
summary(adf_temp) # test-statistics are -0.5237, 0.6261 
# The critical values for the ADF test with a drift are approximately:
# Critical values for test statistics: 
#       1pct  5pct 10pct
# tau2 -3.46 -2.88 -2.57
# phi1  6.52  4.63  3.81
# Both are non-stationary, so we can proceed with cointegration testing.
# Compare the tau2 statistic against the critical values reported in
# the "Critical values for test statistics" table (drift case above).
# Consistent with Schmith et al., both fail to reject the unit-root
# null under this specification.

# ------------------------------------------------------------
# 5. Run FIEG, BIEG, and GIEG on the full sample
# ------------------------------------------------------------
source("R/run_scan_test.R")
COEF_TABLE_FILE <- "output/cv_output/response_surface_coefficients.csv"
R0 <- 0.15

coef_table <- NULL
if (!is.null(COEF_TABLE_FILE)) {
  if (file.exists(COEF_TABLE_FILE)) {
    coef_table <- fread(COEF_TABLE_FILE)
  } else {
    message("coef_table file not found at '", COEF_TABLE_FILE, "' -- running scans without significance.")
  }
}

y <- if (Y_SERIES == "sea_level") combined$sea_level else combined$temperature
X <- if (Y_SERIES == "sea_level") combined$temperature else combined$sea_level

SCHEMES <- c("FIEG", "BIEG", "GIEG")

results <- setNames(
  lapply(SCHEMES, function(scheme) {
    cat(sprintf("\n=== %s (full sample) ===\n", scheme))
    res <- test_cointegration(
      y = y,
      X = X,
      scheme = scheme,
      r0 = R0,
      coef_table = coef_table,
      dates = combined$year,
      min_window = 20
    )
    print(res)
    if (scheme != "GIEG") plot(res)
    res
  }),
  SCHEMES
)

# Access individually, e.g. results$FIEG, results$BIEG, results$GIEG

# ------------------------------------------------------------
# Recursive break search: keep splitting on the GIEG min_t window
# until a resulting subsample is too small to test further.
# ------------------------------------------------------------

min_required <- 10 / R0

run_break_recursive <- function(sub, label = "Full", depth = 0, results_list = list()) {

  indent <- strrep("  ", depth)
  cat(sprintf("\n%s========== [%s] %d obs (%d-%d) ==========\n",
              indent, label, nrow(sub), min(sub$year), max(sub$year)))

  if (nrow(sub) < min_required) {
    message(sprintf("%s%s: only %d obs (< %.0f) -- stopping here, no test run.",
                     indent, label, nrow(sub), min_required))
    results_list[[label]] <- list(data = sub, tests = NULL)
    return(results_list)
  }

  y_sub <- if (Y_SERIES == "sea_level") sub$sea_level else sub$temperature
  X_sub <- if (Y_SERIES == "sea_level") sub$temperature else sub$sea_level

  tests <- setNames(
    lapply(SCHEMES, function(scheme) {
      cat(sprintf("\n%s--- %s ---\n", indent, scheme))
      res <- test_cointegration(
        y = y_sub, X = X_sub, scheme = scheme, r0 = R0,
        coef_table = coef_table, dates = sub$year,
        min_window = 20
      )
      print(res)   # prints: statistic, break window (with dates), critical values + reject verdicts
      res
    }),
    SCHEMES
  )

  results_list[[label]] <- list(data = sub, tests = tests)

  gieg_res <- tests$GIEG

  if (is.na(gieg_res$best_start) || is.na(gieg_res$best_end)) {
    message(sprintf("%s%s: no finite GIEG min_t window -- stopping here.", indent, label))
    return(results_list)
  }

  n_total    <- nrow(sub)
  idx_before <- seq_len(gieg_res$best_start - 1)
  idx_after  <- seq(gieg_res$best_end + 1, n_total)

  can_before <- length(idx_before) >= min_required
  can_after  <- length(idx_after)  >= min_required

  if (!can_before && !can_after) {
    message(sprintf("%s%s: both remaining sides too small to test further -- stopping.", indent, label))
    return(results_list)
  }

  if (can_before) {
    results_list <- run_break_recursive(sub[idx_before], paste0(label, ".Before"), depth + 1, results_list)
  } else if (length(idx_before) > 0) {
    message(sprintf("%s%s.Before: only %d obs -- kept as final segment, not tested.", indent, label, length(idx_before)))
    results_list[[paste0(label, ".Before")]] <- list(data = sub[idx_before], tests = NULL)
  }

  if (can_after) {
    results_list <- run_break_recursive(sub[idx_after], paste0(label, ".After"), depth + 1, results_list)
  } else if (length(idx_after) > 0) {
    message(sprintf("%s%s.After: only %d obs -- kept as final segment, not tested.", indent, label, length(idx_after)))
    results_list[[paste0(label, ".After")]] <- list(data = sub[idx_after], tests = NULL)
  }

  results_list
}

all_segments <- run_break_recursive(combined)

# ------------------------------------------------------------
# Extract all confirmed cointegrated windows from the recursive
# break search, then plot them shaded on the original series.
# ------------------------------------------------------------

extract_windows <- function(results_list, scheme = "GIEG", level = "5%") {

  rbindlist(lapply(names(results_list), function(lbl) {

    seg <- results_list[[lbl]]
    if (is.null(seg$tests)) return(NULL)

    res <- seg$tests[[scheme]]
    if (is.null(res) || is.na(res$best_start)) return(NULL)

    # Only keep windows that were actually significant, if critical
    # values were attached; otherwise keep every detected window.
    if (!is.null(res$reject)) {
      if (!isTRUE(res$reject[[level]])) return(NULL)
    }

    data.table(
      segment    = lbl,
      year_start = res$date_start,
      year_end   = res$date_end,
      statistic  = res$statistic
    )
  }))
}

SIG_LEVEL <- "5%"

windows_dt <- extract_windows(all_segments, scheme = "GIEG", level = SIG_LEVEL)
print(windows_dt)

p_windows <- ggplot() +
  { if (nrow(windows_dt) > 0)
      geom_rect(
        data = windows_dt,
        aes(xmin = year_start, xmax = year_end, ymin = -Inf, ymax = Inf),
        fill = "grey50", alpha = 0.25, inherit.aes = FALSE
      )
  } +
  geom_line(data = combined, aes(x = year, y = sea_level, colour = "Sea level (CSIRO GMSL)")) +
  geom_line(data = combined, aes(x = year, y = rescale_lin(temperature, range_temp, range_sl),
                                  colour = "Temperature anomaly (HadCRUT5)")) +
  { if (nrow(windows_dt) > 0)
      geom_vline(
        data = data.table(edge = c(windows_dt$year_start, windows_dt$year_end)),
        aes(xintercept = edge),
        linetype = "dashed", colour = "black", linewidth = 0.4
      )
  } +
  scale_y_continuous(
    name = "Sea level (mm)",
    sec.axis = sec_axis(
      transform = ~ rescale_lin(., range_sl, range_temp),
      name = "Temperature anomaly (degC, 1850-1900 baseline)"
    )
  ) +
  scale_colour_manual(name = NULL, values = c("steelblue", "firebrick")) +
  labs(
    x = "Year",
    #title = "GIEG-confirmed cointegration windows",
    #subtitle = sprintf("Shaded regions: significant at %s (recursive break search)", SIG_LEVEL)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_windows)

ggsave("output/figures/cointegration_windows.eps",
       p_windows,
       device = cairo_ps,
       width = 8,
       height = 6,
       units = "in")
