# =====================================================================
# paper/gen_surface.R  --  tab_surface.tex (Table 3), standalone
#
#   Rscript paper/gen_surface.R
#
# Environment: OUT [paper/]
# =====================================================================

SIM <- "."
OUT <- Sys.getenv("OUT", "paper")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ndf <- read.csv(file.path(SIM, "results", "tables", "surface_ndf.csv"), stringsAsFactors = FALSE)
TESTS <- c("FIEG", "BIEG", "GIEG"); LV <- c(0.01, 0.05, 0.10)

comma <- function(s) gsub("(?<=\\d)(?=(\\d{3})+($|\\.))", "{,}", s, perl = TRUE)
fmt <- function(x, d) if (is.finite(x)) sprintf("$%s$", comma(sprintf(paste0("%.", d, "f"), x))) else ""
se_ <- function(x, d) if (is.finite(x)) sprintf("(%s)", comma(sprintf(paste0("%.", d, "f"), x))) else ""

main <- ndf[round(ndf$prob, 5) %in% LV, ]
lines <- c()
for (case in c("c", "ct")) {
  panel <- if (case == "c") "A" else "B"
  lines <- c(lines, "\\midrule",
    sprintf("\\multicolumn{8}{l}{\\textit{Panel %s: Case $\\mathrm{%s}$}} \\\\", panel, case), "\\midrule")
  for (tst in TESTS) for (N in 1:3) for (a in LV) {
    r <- main[main$case == case & main$test == tst & main$N == N & abs(main$prob - a) < 1e-9, ]
    if (nrow(r) == 0) next
    r <- r[1, ]
    lbl <- if (N == 1 && a == 0.01) sprintf("\\textit{%s}", tst) else ""
    b3 <- if (r$cubic) fmt(r$b3, 0) else ""; s3 <- if (r$cubic) se_(r$se3, 0) else ""
    lines <- c(lines,
      sprintf("%s & $%d$ & $%d\\%%$ & %s & %s & %s & %s & $%.2f$ \\\\",
              lbl, N, round(100 * a), fmt(r$beta_inf, 3), fmt(r$b1, 2), fmt(r$b2, 1), b3, r$sigma),
      sprintf(" & & & %s & %s & %s & %s & \\\\", se_(r$se_inf, 3), se_(r$se1, 2), se_(r$se2, 1), s3))
  }
}
writeLines(lines, file.path(OUT, "tab_surface.tex"))
cat("ok\n")
