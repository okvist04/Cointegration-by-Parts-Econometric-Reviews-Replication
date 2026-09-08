library(ggplot2)
library(patchwork)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)][1]))
sp   <- Sys.getenv("OUT", here)
gdir <- Sys.getenv("GRAPHICS", file.path(here, "..", "results", "graphics"))
dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

cur <- read.csv(file.path(sp, "fig_ndf_curve.csv"))
pv  <- read.csv(file.path(sp, "fig_ndf_pvals.csv"))
ap  <- read.csv(file.path(sp, "fig_ndf_app.csv"))

lev  <- c("FIEG", "BIEG", "GIEG")
cols <- c(FIEG = "#1b6ca8", BIEG = "#b02418", GIEG = "#2e7d32")
sig  <- c(0.01, 0.05, 0.10)

cur$test <- factor(cur$test, levels = lev)
pv$test  <- factor(pv$test,  levels = lev)
ap$test  <- factor(ap$test,  levels = lev)

# ---------------------------------------------------------------------
# panel (a): fitted distribution function at T = 140, with shaded
# significance bands (darkest = 1%, lightest = 10%) and the
# application's statistics marked as diamonds
# ---------------------------------------------------------------------

bands <- data.frame(
  ymin = c(0, sig[1], sig[2]),
  ymax = c(sig[1], sig[2], sig[3]),
  lvl  = factor(c("1%", "5%", "10%"), levels = c("1%", "5%", "10%"))
)

p_a <- ggplot() +
  geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = lvl),
            inherit.aes = FALSE) +
  scale_fill_manual(name = NULL, values = c("1%" = "grey80", "5%" = "grey87", "10%" = "grey93")) +
  geom_hline(yintercept = sig, colour = "grey45", linetype = "dashed", linewidth = 0.4) +
  geom_line(data = cur, aes(x = stat, y = p, colour = test), linewidth = 0.9) +
  geom_point(data = ap, aes(x = stat, y = pvalue, colour = test), shape = 18, size = 4) +
  scale_colour_manual(name = NULL, values = cols) +
  coord_cartesian(xlim = c(-9, -4), ylim = c(0, 0.25)) +
  labs(x = "Test statistic", y = "p-value") +
  guides(fill = "none") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 10))

# ---------------------------------------------------------------------
# panel (b): calibration on independent null draws -- empirical CDF of
# the assigned p-values against the 45-degree line
# ---------------------------------------------------------------------

p_b <- ggplot(pv, aes(x = p, colour = test)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.5) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  scale_colour_manual(name = NULL, values = cols) +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2)) +
  labs(x = "Nominal p-value", y = "Empirical rejection frequency") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 10))

p_combined <- p_a + p_b + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(gdir, "ndf_dist.eps"), p_a, device = cairo_ps,
       width = 8, height = 6, units = "in")
cat("wrote", file.path(gdir, "ndf.eps"), "\n")

ggsave(file.path(gdir, "ndf_cal.eps"), p_b, device = cairo_ps,
       width = 8, height = 6, units = "in")
cat("wrote", file.path(gdir, "ndf.eps"), "\n")
