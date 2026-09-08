args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)][1]))
sp   <- Sys.getenv("OUT", here)
# Default changed from the original author's sibling Overleaf checkout
# to a folder inside this repository; override with GRAPHICS= to point
# at your actual paper's graphics/ folder.
gdir <- Sys.getenv("GRAPHICS", file.path(here, "..", "results", "graphics"))
dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
cur <- read.csv(file.path(sp,"fig_ndf_curve.csv"))
pv  <- read.csv(file.path(sp,"fig_ndf_pvals.csv"))
ap  <- read.csv(file.path(sp,"fig_ndf_app.csv"))
lev  <- c("FIEG","BIEG","GIEG")
# palette shared with Figure 4 (cointegration_windows.eps)
cols <- c(FIEG="#1b6ca8", BIEG="#b02418", GIEG="#2e7d32")

postscript(file.path(gdir,"ndf.eps"), width=8.8, height=4.4,
           onefile=FALSE, horizontal=FALSE, paper="special", family="Helvetica")
par(mfrow=c(1,2), mar=c(4.2,4.4,2.0,1.2), mgp=c(2.6,0.7,0),
    cex.axis=0.9, cex.lab=1.0, las=1, bty="n")

## panel (a): the fitted distribution function at T = 140
yt <- seq(0, 0.25, by=0.05)
plot(NA, xlim=c(-9,-4), ylim=c(0,0.25), xlab="Test statistic",
     ylab=expression(italic(p)*"-value"), axes=FALSE)
sig <- c(0.01, 0.05, 0.10)
abline(h=setdiff(yt, sig), col="grey92", lwd=0.6)   # plain grid where no level sits
abline(h=sig, col="grey70", lty=2, lwd=0.9)         # the 1%, 5%, 10% levels
axis(1); axis(2, at=yt, labels=format(yt, nsmall=2, trim=TRUE))
for (t in lev) { d <- cur[cur$test==t,]; lines(d$stat, d$p, col=cols[t], lwd=1.8) }
for (t in lev) { a <- ap[ap$test==t,]; points(a$stat, a$pvalue, col=cols[t], pch=18, cex=1.5) }
legend("topleft", legend=lev, col=cols[lev], lwd=1.8, bty="n", cex=0.85, seg.len=1.6)
mtext("(a) Numerical distribution function, Case ct, N = 1, T = 140",
      side=3, line=0.5, cex=0.9)

## panel (b): calibration on 3,000 independent null draws at T = 140
gt <- seq(0, 0.20, by=0.05)
plot(NA, xlim=c(0,0.2), ylim=c(0,0.2), xlab=expression("Nominal "*italic(p)*"-value"),
     ylab="Empirical rejection frequency", axes=FALSE)
abline(h=gt, col="grey92", lwd=0.6)
axis(1, at=gt, labels=format(gt, nsmall=2, trim=TRUE))
axis(2, at=gt, labels=format(gt, nsmall=2, trim=TRUE))
abline(0, 1, col="grey40", lty=2, lwd=1.1)
for (t in lev) {
  x <- sort(pv$p[pv$test==t]); F <- seq_along(x)/length(x)
  k <- x <= 0.205
  lines(x[k], F[k], col=cols[t], lwd=1.8, type="s")
}
legend("topleft", legend=lev, col=cols[lev], lwd=1.8, bty="n", cex=0.85, seg.len=1.6)
mtext("(b) Calibration on 3,000 independent null draws", side=3, line=0.5, cex=0.9)
dev.off()
cat("wrote", file.path(gdir,"ndf.eps"), "\n")
