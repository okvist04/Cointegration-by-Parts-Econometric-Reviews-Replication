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
d <- merge(sl[, c("year", "GMSL (mm)")], tp[, c("year", "anom")], by = "year")
d <- d[d$year >= YEAR_MIN & d$year <= YEAR_MAX, ]
names(d)[2] <- "gmsl"
d <- d[order(d$year), ]

colS <- "#1b6ca8"; colT <- "#b02418"

postscript(file.path(G, "GMST_GMTA.eps"), width = 8.4, height = 4.8,
           onefile = FALSE, horizontal = FALSE, paper = "special", family = "Helvetica")
par(mar = c(4.2, 4.4, 1.2, 4.6), mgp = c(2.6, 0.7, 0), cex.axis = 0.9, las = 1, bty = "n")
plot(NA, xlim = range(d$year), ylim = range(d$gmsl), xlab = "Year", ylab = "Sea level (mm)", axes = FALSE)
abline(h = pretty(d$gmsl), col = "grey92", lwd = 0.6)
axis(1); axis(2)
lines(d$year, d$gmsl, col = colS, lwd = 1.6)

# right axis for temperature (linear rescale onto sea level's y-range
# purely for display -- the vertical alignment matches RANGES, not any
# principled unit conversion, so visual closeness between the two lines
# isn't itself meaningful)
r <- range(d$anom); u <- par("usr")[3:4]
sc <- function(z) u[1] + (z - r[1]) / diff(r) * diff(u)
lines(d$year, sc(d$anom), col = colT, lwd = 1.6)
at <- pretty(r); at <- at[at >= r[1] & at <= r[2]]
axis(4, at = sc(at), labels = format(at, trim = TRUE))
mtext(expression("Temperature anomaly (" * degree * "C, 1850-1900 baseline)"), side = 4, line = 3, las = 0, cex = 1.0)
legend("topleft", legend = c("Sea level (CSIRO GMSL)", "Temperature anomaly (HadCRUT5)"),
       col = c(colS, colT), lwd = 1.6, bty = "n", cex = 0.85)
dev.off()
cat("wrote GMST_GMTA.eps,", YEAR_MIN, "-", YEAR_MAX, "\n")
