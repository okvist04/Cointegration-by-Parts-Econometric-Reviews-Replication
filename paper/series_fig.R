# =====================================================================
# paper/series_fig.R  --  GMST_GMTA.eps: raw GMSL and temperature
# anomaly series on a dual axis, over the full application sample.
#
#   Rscript paper/series_fig.R
#
# Base R only (postscript device), matching windows_fig.R's style --
# no ggplot2/cairo_ps dependency. Reads data/ directly, same rebasing
# (HadCRUT5 to a 1850-1900 baseline) as calibration/run_application.R.
#
# Environment:
#   APPLICATION_DATA [../data]   GRAPHICS [../results/graphics]
#   YEAR_MIN [1880]  YEAR_MAX [2019]
# =====================================================================
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

# linear rescale of temperature onto sea level's range, purely for
# display on a shared panel -- the right-hand axis back-transforms to
# real units, so the apparent closeness of the two lines isn't itself
# meaningful (a visualization choice, not a result)
range_sl   <- range(combined$sea_level, na.rm = TRUE)
range_temp <- range(combined$temperature, na.rm = TRUE)
rescale_lin <- function(x, from, to) (x - from[1]) / diff(from) * diff(to) + to[1]
combined$temperature_scaled <- rescale_lin(combined$temperature, range_temp, range_sl)

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
    subtitle = sprintf("%d-%d, temperature rescaled onto sea level's range for display",
                        YEAR_MIN, YEAR_MAX)
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(G, "GMST_GMTA.eps"), p_overlay, device = cairo_ps,
       width = 8, height = 6, units = "in")
cat("wrote", file.path(G, "GMST_GMTA.eps"), "\n")
