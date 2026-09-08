# =====================================================================
# paper/ndf_fig_data.R  --  fig_ndf_*.csv, consumed by ndf_fig.R
#
#   Rscript paper/ndf_fig_data.R
#
# Environment: OUT [paper/]
# =====================================================================

source("R/response_surface.R")

SIM <- "."
OUT <- Sys.getenv("OUT", "paper")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ndf <- read.csv(file.path(SIM, "results", "tables", "surface_ndf.csv"), stringsAsFactors = FALSE)
TESTS <- c("FIEG", "BIEG", "GIEG"); CASE <- "ct"; T <- 140; N <- 1

sfor <- function(tst) {
  s <- ndf[ndf$case == CASE & ndf$N == N & ndf$test == tst, ]
  s <- s[order(s$prob), ]
  list(probs = s$prob, surfaces = lapply(seq_len(nrow(s)), function(i) {
    r <- s[i, ]
    list(beta = c(r$beta_inf, r$b1, r$b2, r$b3), se = c(r$se_inf, r$se1, r$se2, r$se3),
         sigma = r$sigma, cubic = r$cubic, dof = r$dof, nobs = r$nobs)
  }))
}

# panel (a): fitted distribution function at T = 140
curve_rows <- list()
for (tst in TESTS) {
  sf <- sfor(tst)
  cvals <- vapply(sf$surfaces, cv_at, numeric(1), T = T)
  xs <- seq(min(cvals), -2.0, length.out = 600)
  for (x in xs) {
    curve_rows[[length(curve_rows) + 1]] <- data.frame(test = tst, stat = x,
      p = pvalue_from_surface(sf$probs, sf$surfaces, T, x))
  }
}
write.csv(do.call(rbind, curve_rows), file.path(OUT, "fig_ndf_curve.csv"), row.names = FALSE)

# panel (b): p-values of independent null draws at T = 140
st <- readRDS(file.path(SIM, "results", "sizepower", sprintf("size_case-%s_T-%d_N-%d.rds", CASE, T, N)))$stats
pv_rows <- list()
for (tst in TESTS) {
  sf <- sfor(tst)
  v <- st[[tolower(tst)]]
  for (x in v) {
    pv_rows[[length(pv_rows) + 1]] <- data.frame(test = tst, p = pvalue_from_surface(sf$probs, sf$surfaces, T, x))
  }
}
write.csv(do.call(rbind, pv_rows), file.path(OUT, "fig_ndf_pvals.csv"), row.names = FALSE)

# application statistics
app <- read.csv(file.path(SIM, "results", "application", "application.csv"), stringsAsFactors = FALSE)
h <- app[app$case == CASE & app$r0 == 0.15 & app$criterion == "BIC" & app$normalization == "temp on sl",
         c("test", "stat", "pvalue")]
write.csv(h, file.path(OUT, "fig_ndf_app.csv"), row.names = FALSE)
cat("ok\n")
