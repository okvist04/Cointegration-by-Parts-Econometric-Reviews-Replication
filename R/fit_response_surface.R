# ============================================================
# fit_response_surface.R
#
# Turns the simulated critical values into a MacKinnon (2010)-style
# response surface: for each (scan_type, n_regressors, prob), fit
#
#   crit_val = b0 + b1/T + b2/T^2 [+ b3/T^3]
#
# by weighted least squares, weighting each simulated T by the inverse
# variance of its own critical-value estimate (1 / crit_val_se^2) --
# this is the standard MacKinnon weighting, and it's exactly what the
# bootstrap SE already saved per combination is for.
#
# THREE ways to get a coefficient table, depending on how you want to
# decide between quadratic (max_power=2) and cubic (max_power=3):
#   - fit_response_surfaces(results, max_power):    one power for every
#     group, chosen by you up front.
#   - diagnose_response_surface_power(results):     doesn't fit a final
#     table -- fits both powers per group and reports whether the cubic
#     term is statistically significant, so you can decide with actual
#     evidence rather than a blanket assumption. Run this FIRST.
#   - fit_response_surfaces_auto(results, alpha):   fits both per group
#     and automatically keeps the cubic term only where it's
#     significant at the given alpha -- MacKinnon's own general-to-
#     specific approach (1994, 1996), applied per group.
#
# Output:
#   - response_surface_coefficients.csv: one row per
#     (scan_type, n_regressors, prob), with b0, b1, b2, b3 (NA if
#     unused), residual info, and the T-range the fit was built from.
#   - get_critical_value(): looks up the right row and evaluates the
#     surface at any T (see examples at the bottom of this file).
#
# Run this AFTER run_generation_cv.R has produced
# critical_values_null_full_grid.csv (or after enough checkpoints
# exist in checkpoints/ that you want to fit against what's done so far).
# ============================================================

library(data.table)

# ------------------------------------------------------------
# Internal helpers (shared by all three fitting functions below, so
# the actual regression logic -- weights, formula, which columns end
# up in a result row -- exists in exactly one place).
# ------------------------------------------------------------

.compute_weights <- function(sub) {
  # Weight by inverse variance of the empirical quantile estimate.
  # Guard against zero/NA SE (shouldn't happen, but a WLS weight of
  # Inf would silently break the fit rather than error).
  se <- sub$crit_val_se
  se[!is.finite(se) | se <= 0] <- NA_real_
  if (all(is.na(se))) {
    return(rep(1, nrow(sub)))
  }
  se[is.na(se)] <- max(se, na.rm = TRUE)  # fall back to least-informative weight
  1 / se^2
}

.fit_power <- function(sub, power) {
  w <- .compute_weights(sub)
  inv_T <- 1 / sub$T
  inv_T2 <- inv_T^2
  if (power == 2) {
    lm(crit_val ~ inv_T + inv_T2, data = sub, weights = w)
  } else {
    inv_T3 <- inv_T^3
    lm(crit_val ~ inv_T + inv_T2 + inv_T3, data = sub, weights = w)
  }
}

.surface_row_from_fit <- function(fit, st_i, k_i, p_i, sub) {
  cf <- coef(fit)
  s <- summary(fit)
  b3 <- if ("inv_T3" %in% names(cf)) cf[["inv_T3"]] else NA_real_
  data.table(
    scan_type = st_i,
    n_regressors = k_i,
    prob = p_i,
    b0 = cf[["(Intercept)"]],
    b1 = cf[["inv_T"]],
    b2 = cf[["inv_T2"]],
    b3 = b3,
    resid_se = s$sigma,
    r_squared = s$r.squared,
    n_T_points = nrow(sub),
    T_min = min(sub$T),
    T_max = max(sub$T)
  )
}

.usable_groups <- function(results) {
  results <- as.data.table(results)
  list(
    results = results,
    groups = unique(results[, .(scan_type, n_regressors, prob)])
  )
}

# ------------------------------------------------------------
# fit_response_surfaces
# ------------------------------------------------------------
#
# results: data.table with columns T, n_regressors, scan_type, prob,
#          crit_val, crit_val_se (as produced by
#          calibrate_null_critical_values_parallel()).
# max_power: highest power of 1/T to include (2 or 3), applied to
#          EVERY group. Default 2 is the safer choice unless you have
#          many T points per group -- with 12 T values per
#          (scan_type, n_regressors, prob) group, a cubic term (4
#          parameters) is usable but leaves few residual degrees of
#          freedom. Use diagnose_response_surface_power() first if you
#          want evidence for this choice rather than an assumption, or
#          fit_response_surfaces_auto() to decide per group instead of
#          picking one power for the whole table.

fit_response_surfaces <- function(
    results,
    max_power = 2
) {
  
  if (!max_power %in% c(2, 3))
    stop("max_power must be 2 or 3.")
  
  u <- .usable_groups(results)
  results <- u$results
  groups <- u$groups
  
  out <- vector("list", nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    
    st_i <- groups$scan_type[i]
    k_i <- groups$n_regressors[i]
    p_i <- groups$prob[i]
    
    sub <- results[
      scan_type == st_i & n_regressors == k_i & prob == p_i
    ][order(T)]
    
    sub <- sub[is.finite(crit_val)]
    
    if (nrow(sub) < max_power + 2) {
      warning(sprintf(
        "Skipping scan_type=%s n_regressors=%d prob=%s: only %d usable T points.",
        st_i, k_i, p_i, nrow(sub)
      ))
      next
    }
    
    fit <- .fit_power(sub, max_power)
    out[[i]] <- .surface_row_from_fit(fit, st_i, k_i, p_i, sub)
  }
  
  rbindlist(out)
}

# ------------------------------------------------------------
# diagnose_response_surface_power
# ------------------------------------------------------------
#
# Run this BEFORE choosing max_power. For each (scan_type,
# n_regressors, prob) group, fits BOTH the quadratic and cubic
# surfaces and reports whether the cubic term (b3) is statistically
# significant -- i.e. whether the data actually support the extra
# curvature, rather than just adding a parameter that fits noise.
#
# This does not produce a usable coef_table on its own -- it's a
# diagnostic to inform the max_power choice (or to report directly in
# the paper as justification for whichever choice you make).
#
# alpha: significance threshold for flagging the cubic term as needed
#        (default 0.05, i.e. |t| large enough for p < 0.05).
#
# Returns one row per group with:
#   b3_estimate, b3_se, b3_tstat, b3_pvalue  -- the test itself
#   cubic_significant                        -- TRUE/FALSE at `alpha`
#   resid_se_quadratic, resid_se_cubic       -- compare fit quality
#   r_squared_quadratic, r_squared_cubic
#   n_T_points

diagnose_response_surface_power <- function(
    results,
    alpha = 0.05
) {
  
  u <- .usable_groups(results)
  results <- u$results
  groups <- u$groups
  
  out <- vector("list", nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    
    st_i <- groups$scan_type[i]
    k_i <- groups$n_regressors[i]
    p_i <- groups$prob[i]
    
    sub <- results[
      scan_type == st_i & n_regressors == k_i & prob == p_i
    ][order(T)]
    
    sub <- sub[is.finite(crit_val)]
    
    # Need at least 5 points to fit a cubic (4 parameters) with any
    # residual degrees of freedom at all.
    if (nrow(sub) < 5) {
      warning(sprintf(
        "Skipping scan_type=%s n_regressors=%d prob=%s: only %d usable T points (need >= 5 to test the cubic term).",
        st_i, k_i, p_i, nrow(sub)
      ))
      next
    }
    
    fit_quad <- .fit_power(sub, 2)
    fit_cubic <- .fit_power(sub, 3)
    
    s_quad <- summary(fit_quad)
    s_cubic <- summary(fit_cubic)
    
    cf_cubic <- s_cubic$coefficients
    b3_est <- cf_cubic["inv_T3", "Estimate"]
    b3_se <- cf_cubic["inv_T3", "Std. Error"]
    b3_t <- cf_cubic["inv_T3", "t value"]
    b3_p <- cf_cubic["inv_T3", "Pr(>|t|)"]
    
    out[[i]] <- data.table(
      scan_type = st_i,
      n_regressors = k_i,
      prob = p_i,
      n_T_points = nrow(sub),
      b3_estimate = b3_est,
      b3_se = b3_se,
      b3_tstat = b3_t,
      b3_pvalue = b3_p,
      cubic_significant = is.finite(b3_p) & b3_p < alpha,
      resid_se_quadratic = s_quad$sigma,
      resid_se_cubic = s_cubic$sigma,
      r_squared_quadratic = s_quad$r.squared,
      r_squared_cubic = s_cubic$r.squared
    )
  }
  
  rbindlist(out)
}

# ------------------------------------------------------------
# fit_response_surfaces_auto
# ------------------------------------------------------------
#
# Like fit_response_surfaces(), but decides quadratic vs. cubic PER
# GROUP rather than applying one power to the whole table: fits the
# cubic term, keeps it only if significant at `alpha`, otherwise falls
# back to quadratic. This is MacKinnon's own general-to-specific
# approach (1994, 1996) -- fit the richer specification, drop terms
# that aren't significant -- applied separately to each
# (scan_type, n_regressors, prob) group rather than the table as a
# whole.
#
# Adds a `power_used` column (2 or 3) to the output so you can see,
# and report, which groups needed the cubic term -- e.g. summarizing
# with table(coef_table$power_used, coef_table$scan_type) to check for
# a pattern (such as GIEG needing it more often than FIEG/BIEG, which
# would itself be worth a sentence in the paper).

fit_response_surfaces_auto <- function(
    results,
    alpha = 0.05
) {
  
  u <- .usable_groups(results)
  results <- u$results
  groups <- u$groups
  
  out <- vector("list", nrow(groups))
  
  for (i in seq_len(nrow(groups))) {
    
    st_i <- groups$scan_type[i]
    k_i <- groups$n_regressors[i]
    p_i <- groups$prob[i]
    
    sub <- results[
      scan_type == st_i & n_regressors == k_i & prob == p_i
    ][order(T)]
    
    sub <- sub[is.finite(crit_val)]
    
    if (nrow(sub) < 4) {
      warning(sprintf(
        "Skipping scan_type=%s n_regressors=%d prob=%s: only %d usable T points.",
        st_i, k_i, p_i, nrow(sub)
      ))
      next
    }
    
    if (nrow(sub) >= 5) {
      fit_cubic <- .fit_power(sub, 3)
      b3_p <- summary(fit_cubic)$coefficients["inv_T3", "Pr(>|t|)"]
      use_cubic <- is.finite(b3_p) && b3_p < alpha
    } else {
      # Not enough points to even fit the cubic term -- quadratic only.
      use_cubic <- FALSE
    }
    
    power_used <- if (use_cubic) 3L else 2L
    fit <- .fit_power(sub, power_used)
    row <- .surface_row_from_fit(fit, st_i, k_i, p_i, sub)
    row[, power_used := power_used]
    out[[i]] <- row
  }
  
  rbindlist(out)
}

# ------------------------------------------------------------
# get_critical_value
# ------------------------------------------------------------
#
# Looks up the fitted surface for the requested (scan_type,
# n_regressors, prob) and evaluates it at the requested sample size T.
# Warns (does not error) if T falls outside the range of simulated T
# values the surface was fit on, since the polynomial-in-1/T form is
# only validated to interpolate well, not extrapolate.

get_critical_value <- function(
    T,
    n_regressors,
    scan_type,
    prob,
    coef_table
) {

  match_idx <- coef_table$scan_type == scan_type &
    coef_table$n_regressors == n_regressors &
    coef_table$prob == prob

  row <- coef_table[match_idx]

  if (nrow(row) == 0) {
    stop(sprintf(
      "No fitted surface for scan_type=%s, n_regressors=%s, prob=%s.",
      scan_type, n_regressors, prob
    ))
  }
  if (nrow(row) > 1) {
    stop("Multiple matching rows in coef_table -- coef_table should have one row per (scan_type, n_regressors, prob).")
  }

  if (T < row$T_min || T > row$T_max) {
    warning(sprintf(
      "T = %d is outside the simulated range [%d, %d] for this combination -- extrapolating.",
      T, row$T_min, row$T_max
    ))
  }

  val <- row$b0 + row$b1 / T + row$b2 / T^2
  if (!is.na(row$b3)) {
    val <- val + row$b3 / T^3
  }

  val
}

# ============================================================
# Example usage (uncomment to run after the full grid is done):
#
OUTPUT_DIR <- "output/cv_output"   # must match the folder used in
#                              # run_generation_cv.R / validate_fast_vs_r.R
#
results <- fread(file.path(OUTPUT_DIR, "critical_values_null_full_grid.csv"))
#
# --- Step 1: decide on max_power with actual evidence ---
diagnostics <- diagnose_response_surface_power(results, alpha = 0.05)
print(diagnostics[order(-cubic_significant, scan_type, n_regressors, prob)])
# # How many groups actually need the cubic term?
table(diagnostics$cubic_significant, diagnostics$scan_type)
#
# --- Step 2a: EITHER pick one power for the whole table ... ---
# coef_table <- fit_response_surfaces(results, max_power = 2)
#
# --- Step 2b: ... OR let it decide per group automatically ---
coef_table <- fit_response_surfaces_auto(results, alpha = 0.05)
table(coef_table$power_used, coef_table$scan_type)  # which groups used cubic
#
fwrite(coef_table, file.path(OUTPUT_DIR, "response_surface_coefficients.csv"))
#
# get_critical_value(
#   T = 623, n_regressors = 3, scan_type = "all", prob = 0.05,
#   coef_table = coef_table
# )
# ============================================================
