library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)][1]))
D <- Sys.getenv("APPLICATION_DATA", file.path(here, "..", "data"))
G <- Sys.getenv("GRAPHICS", file.path(here, "..", "results", "graphics"))
dir.create(G, recursive = TRUE, showWarnings = FALSE)
YEAR_MIN <- as.integer(Sys.getenv("YEAR_MIN", "1880"))
YEAR_MAX <- as.integer(Sys.getenv("YEAR_MAX", "2019"))

sl <- read.csv(file.path(D, "CSIRO_Recons_gmsl_yr_2019.csv"), check.names = FALSE)
tp <- read.csv(file.path(D, "HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"), check.names = FALSE)
sl$year <- floor(sl$Time); tp$year <- floor(tp$Time)
base <- mean(tp[["Anomaly (deg C)"]][tp$year >= 1850 & tp$year <= 1900])
tp$anom <- tp[["Anomaly (deg C)"]] - base
combined <- merge(sl[, c("year", "GMSL (mm)")], tp[, c("year", "anom")], by = "year")
names(combined) <- c("year", "sea_level", "temperature")
combined <- combined[combined$year >= YEAR_MIN & combined$year <= YEAR_MAX, ]
combined <- combined[order(combined$year), ]

range_sl   <- range(combined$sea_level, na.rm = TRUE)
range_temp <- range(combined$temperature, na.rm = TRUE)
rescale_lin <- function(x, from, to) (x - from[1]) / diff(from) * diff(to) + to[1]
combined$temperature_scaled <- rescale_lin(combined$temperature, range_temp, range_sl)

# Window read from the actual run_application.R output (headline
# specification: r0 = 0.15, BIC, temperature on sea level, case ct) --
# override via env vars for a different row's window.
app_csv <- Sys.getenv("APPLICATION_CSV", "results/application/application.csv")
if (!file.exists(app_csv)) {
  stop(sprintf(
    "%s not found -- run calibration/run_application.R first (it writes the GIEG window this figure shades).",
    app_csv))
}
app <- read.csv(app_csv, stringsAsFactors = FALSE)
WIN_CASE <- Sys.getenv("WIN_CASE", "ct")
WIN_R0 <- as.numeric(Sys.getenv("WIN_R0", "0.15"))
WIN_CRITERION <- Sys.getenv("WIN_CRITERION", "BIC")
WIN_NORM <- Sys.getenv("WIN_NORM", "temp on sl")
row <- app[app$case == WIN_CASE & abs(app$r0 - WIN_R0) < 1e-9 &
           app$criterion == WIN_CRITERION & app$normalization == WIN_NORM &
           app$test == "GIEG", ]
if (nrow(row) != 1) {
  stop(sprintf("expected exactly one matching GIEG row in %s for case=%s r0=%s criterion=%s normalization=%s, found %d",
               app_csv, WIN_CASE, WIN_R0, WIN_CRITERION, WIN_NORM, nrow(row)))
}
win_start <- row$win_start[1]; win_end <- row$win_end[1]

p_windows <- ggplot() +
  geom_rect(aes(xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf),
            fill = "grey50", alpha = 0.25) +
  geom_line(data = combined, aes(x = year, y = sea_level, colour = "Sea level (CSIRO GMSL)")) +
  geom_line(data = combined, aes(x = year, y = temperature_scaled, colour = "Temperature anomaly (HadCRUT5)")) +
  geom_vline(xintercept = c(win_start, win_end), linetype = "dashed", colour = "black", linewidth = 0.4) +
  scale_y_continuous(
    name = "Sea level (mm)",
    sec.axis = sec_axis(
      transform = ~ rescale_lin(., range_sl, range_temp),
      name = "Temperature anomaly (degC, 1850-1900 baseline)"
    )
  ) +
  scale_colour_manual(name = NULL, values = c("steelblue", "firebrick")) +
  labs(x = "Year") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(G, "cointegration_windows.eps"), p_windows, device = cairo_ps,
       width = 8, height = 6, units = "in")
cat("wrote cointegration_windows.eps, window", win_start, "-", win_end, "\n")
