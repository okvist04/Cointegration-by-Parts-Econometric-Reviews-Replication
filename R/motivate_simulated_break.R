# ============================================================
# motivate_simulated_break.R
#
# R equivalent of Example_simulated_data.qmd: simulates data with two
# regime breaks (coint -> noncoint -> coint), then shows that GIEG(r0)'s
# scan across subsamples correctly picks out the cointegrated segment
# as its minimum-t window, despite the breaks elsewhere in the sample.
#
# Uses simulate_window_break() (generate_cv.R) and run_scan_test()
# (run_scan_test.R, return_full_scan = TRUE) -- no new machinery needed.
#
# Needs the "patchwork" package to combine the two panels side by side
# (install.packages("patchwork") if you don't have it) -- ggplot2 alone
# doesn't arrange separate plot objects next to each other.
# ============================================================

source("R/generate_cv.R")
source("R/run_scan_test.R")
library(patchwork)

# ------------------------------------------------------------
# 1. Simulate data with two breaks
# ------------------------------------------------------------
# Matches the co-author's Julia parameters exactly: noncointegrated for
# the first quarter and last quarter of the sample, cointegrated in the
# middle half.

T           <- 200
break_fracs <- c(0.25, 0.75)   # noncoint [1,50] -> coint [51,150] -> noncoint [151,200]
seed        <- 21

dgp <- simulate_window_break(
  T               = T,
  break_fracs     = break_fracs,
  alpha           = 1.5,
  beta            = 1.2,
  sigma_x         = 1.0,
  sigma_e         = 1.5,
  rho_stationary  = 0.4,
  starting_regime = "noncoint",
  seed            = seed
)

break_points <- floor(T * break_fracs)

# ------------------------------------------------------------
# 2. Plot the simulated data
# ------------------------------------------------------------

data_df <- data.frame(time = 1:T, X = as.numeric(dgp$X), Y = dgp$y)

p_data <- ggplot(data_df, aes(x = time)) +
  geom_line(aes(y = X, colour = "X")) +
  geom_line(aes(y = Y, colour = "Y")) +
  geom_vline(xintercept = break_points, linetype = "dashed", colour = "black") +
  scale_colour_manual(name = NULL, values = c(X = "purple", Y = "darkgreen")) +
  labs(x = "Time", y = "Value", title = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 3. Run the GIEG(r0) scan
# ------------------------------------------------------------
# min_window_obs = 40 at T = 200 -> r0 = 0.2, matching the co-author's
# min_obs = min_window_obs = 40 in the Julia version exactly.

full <- run_scan_test(
  y = dgp$y,
  X = dgp$X,
  scheme = "GIEG",
  r0 = 0.15,
  include_intercept = TRUE,
  max_lag = NULL,
  return_full_scan = TRUE
)

cat(sprintf("  Min t-stat : %.4f\n", full$statistic))
cat(sprintf(
  "  Window     : [%d, %d]  (n=%d)\n",
  full$best_start, full$best_end, full$best_end - full$best_start + 1
))

# ------------------------------------------------------------
# 4. Plot every scanned window's t-statistic, with the minimum
#    highlighted -- the "spaghetti" plot showing the scan finds the
#    cointegrated segment despite the breaks elsewhere.
# ------------------------------------------------------------

min_row <- data.frame(start = full$best_start, end = full$best_end, tstat = full$statistic)

p_scan <- ggplot() +
  geom_segment(
    data = full$scan,
    aes(x = start, xend = end, y = tstat, yend = tstat),
    colour = "blue", alpha = 0.2, linewidth = 0.3
  ) +
  geom_segment(
    data = min_row,
    aes(x = start, xend = end, y = tstat, yend = tstat, colour = "Min t-stat"),
    linewidth = 1
  ) +
  geom_vline(
    aes(xintercept = break_points, colour = "Breaks"),
    linetype = "dashed"
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Min t-stat" = "red", "Breaks" = "black")
  ) +
  labs(x = "Time", y = "t-statistic", title = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 5. Combine side by side and save
# ------------------------------------------------------------

print(p_data)
print(p_scan)

ggsave("output/figures/sim_breaks.eps",
       p_data,
       device = cairo_ps,
       width = 8,
       height = 6,
       units = "in")

ggsave("output/figures/t_stat_sim_breaks.eps",
       p_scan,
       device = cairo_ps,
       width = 8,
       height = 6,
       units = "in")
