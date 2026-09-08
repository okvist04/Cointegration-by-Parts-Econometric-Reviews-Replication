args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)][1]))
D <- Sys.getenv("APPLICATION_DATA", file.path(here, "..", "data"))
G <- Sys.getenv("GRAPHICS", file.path(here, "..", "results", "graphics"))
dir.create(G, recursive = TRUE, showWarnings = FALSE)
sl <- read.csv(file.path(D,"CSIRO_Recons_gmsl_yr_2019.csv"), check.names=FALSE)
tp <- read.csv(file.path(D,"HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"), check.names=FALSE)
sl$year <- floor(sl$Time); tp$year <- floor(tp$Time)
base <- mean(tp[["Anomaly (deg C)"]][tp$year >= 1850 & tp$year <= 1900])
tp$anom <- tp[["Anomaly (deg C)"]] - base
d <- merge(sl[,c("year","GMSL (mm)")], tp[,c("year","anom")], by="year")
d <- d[d$year >= 1880 & d$year <= 2019,]
names(d)[2] <- "gmsl"

# Window read from the actual run_application.R output (headline
# specification: r0 = 0.15, BIC, temperature on sea level, case ct) --
# no longer a hardcoded value. Override the specification via the
# WIN_CASE / WIN_R0 / WIN_CRITERION / WIN_NORM env vars if you want a
# different row's window shaded instead.
app_csv <- Sys.getenv("APPLICATION_CSV", file.path(here, "..", "results", "application", "application.csv"))
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
win <- c(row$win_start[1], row$win_end[1])
colS <- "#1b6ca8"; colT <- "#b02418"

postscript(file.path(G,"cointegration_windows.eps"), width=8.4, height=4.8,
           onefile=FALSE, horizontal=FALSE, paper="special", family="Helvetica")
par(mar=c(4.2,4.4,1.2,4.6), mgp=c(2.6,0.7,0), cex.axis=0.9, las=1, bty="n")
plot(NA, xlim=range(d$year), ylim=range(d$gmsl), xlab="Year", ylab="Sea level (mm)", axes=FALSE)
rect(win[1], par("usr")[3], win[2], par("usr")[4], col="grey88", border=NA)
abline(v=win, lty=2, lwd=1.2)
abline(h=pretty(d$gmsl), col="grey92", lwd=0.6)
axis(1); axis(2)
lines(d$year, d$gmsl, col=colS, lwd=1.6)
# right axis for temperature
r <- range(d$anom); u <- par("usr")[3:4]
sc <- function(z) u[1] + (z - r[1])/diff(r) * diff(u)
lines(d$year, sc(d$anom), col=colT, lwd=1.6)
at <- pretty(r); at <- at[at >= r[1] & at <= r[2]]
axis(4, at=sc(at), labels=format(at, trim=TRUE))
mtext(expression("Temperature anomaly ("*degree*"C, 1850-1900 baseline)"), side=4, line=3, las=0, cex=1.0)
legend("topleft", legend=c("Sea level (CSIRO GMSL)","Temperature anomaly (HadCRUT5)"),
       col=c(colS,colT), lwd=1.6, bty="n", cex=0.85)
dev.off()
cat("wrote cointegration_windows.eps, window", win[1], "-", win[2], "\n")
