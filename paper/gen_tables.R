# =====================================================================
# paper/gen_tables.R  --  tab_cv_{c,ct}, tab_ndfcal, tab_app, tab_robust
#
#   Rscript paper/gen_tables.R
#
# Run after make_tables.R, run_size_power.R, run_stage2_null.R, and
# run_application.R have written their outputs under results/.
#
# Environment: OUT [paper/]
# =====================================================================

source("R/response_surface.R")

SIM <- "."
OUT <- Sys.getenv("OUT", "paper")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

q   <- read.csv(file.path(SIM, "results", "tables", "quantiles.csv"), stringsAsFactors = FALSE)
ndf <- read.csv(file.path(SIM, "results", "tables", "surface_ndf.csv"), stringsAsFactors = FALSE)
TESTS <- c("FIEG", "BIEG", "GIEG"); LV <- c(0.01, 0.05, 0.10)

Tfmt <- function(T) {
  if (T >= 1000) {
    s <- sprintf("%d", T)
    sub("^(\\d)(\\d{3})$", "\\1{,}\\2", s)
  } else as.character(T)
}

# ---- Table 4.1 / 4.2 : simulated quantiles N = 1, with MCSE -------
for (case in c("c", "ct")) {
  lines <- c()
  sub <- q[q$case == case & q$N == 1 & round(q$prob, 5) %in% LV, ]
  for (T in sort(unique(sub$T))) {
    v <- c(); s <- c()
    for (tst in TESTS) for (a in LV) {
      r <- sub[sub$T == T & sub$test == tst & abs(sub$prob - a) < 1e-9, ]
      v <- c(v, sprintf("$%.3f$", r$cv[1])); s <- c(s, sprintf("(%.3f)", r$mcse[1]))
    }
    lines <- c(lines,
      paste0("    $", Tfmt(T), "$ & ", paste(v, collapse = " & "), " \\\\"),
      paste0("     & ", paste(s, collapse = " & "), " \\\\[2pt]"))
  }
  writeLines(lines, file.path(OUT, sprintf("tab_cv_%s.tex", case)))
}

# ---- Table 4.3 : response surface, both cases, N = 1..3 -----------
main <- ndf[round(ndf$prob, 5) %in% LV, ]
lines <- c()
for (case in c("c", "ct")) {
  panel <- if (case == "c") "A" else "B"
  lines <- c(lines,
    "\\midrule",
    sprintf("\\multicolumn{8}{l}{\\textit{Panel %s: Case $\\mathrm{%s}$}} \\\\", panel, case),
    "\\midrule")
  for (tst in TESTS) for (N in 1:3) for (a in LV) {
    r <- main[main$case == case & main$test == tst & main$N == N & abs(main$prob - a) < 1e-9, ]
    if (nrow(r) == 0) next
    r <- r[1, ]
    comma <- function(s) gsub("(?<=\\d)(?=(\\d{3})+\\.)", "{,}", s, perl = TRUE)
    f <- function(x) if (is.finite(x)) sprintf("$%s$", comma(sprintf("%.3f", x))) else ""
    g <- function(x) if (is.finite(x)) sprintf("(%s)", comma(sprintf("%.3f", x))) else ""
    lbl <- if (N == 1 && a == 0.01) sprintf("\\textit{%s}", tst) else ""
    b3 <- if (r$cubic) f(r$b3) else ""; s3 <- if (r$cubic) g(r$se3) else ""
    lines <- c(lines,
      sprintf("%s & $%d$ & $%d\\%%$ & %s & %s & %s & %s & $%.2f$ \\\\",
              lbl, N, round(100 * a), f(r$beta_inf), f(r$b1), f(r$b2), b3, r$sigma),
      sprintf(" & & & %s & %s & %s & %s & \\\\", g(r$se_inf), g(r$se1), g(r$se2), s3))
  }
}
writeLines(lines, file.path(OUT, "tab_surface.tex"))

surfaces_for <- function(case, N, tst) {
  s <- ndf[ndf$case == case & ndf$N == N & ndf$test == tst, ]
  s <- s[order(s$prob), ]
  list(probs = s$prob, surfaces = lapply(seq_len(nrow(s)), function(i) {
    r <- s[i, ]
    list(beta = c(r$beta_inf, r$b1, r$b2, r$b3), se = c(r$se_inf, r$se1, r$se2, r$se3),
         sigma = r$sigma, cubic = r$cubic, dof = r$dof, nobs = r$nobs)
  }))
}

# ---- Table 4.4 : NDF calibration on independent null draws --------
lines <- c()
for (case in c("c", "ct")) {
  panel <- if (case == "c") "A" else "B"
  lines <- c(lines,
    "\\midrule",
    sprintf("\\multicolumn{11}{l}{\\textit{Panel %s: Case $\\mathrm{%s}$}} \\\\", panel, case),
    "\\midrule")
  for (T in c(90, 140, 350, 600)) for (N in c(1, 3)) {
    f <- file.path(SIM, "results", "sizepower", sprintf("size_case-%s_T-%d_N-%d.rds", case, T, N))
    if (!file.exists(f)) next
    st <- readRDS(f)$stats
    cells <- c()
    for (tst in TESTS) {
      sf <- surfaces_for(case, N, tst)
      v <- st[[tolower(tst)]]
      pv <- vapply(v, function(x) pvalue_from_surface(sf$probs, sf$surfaces, T, x), numeric(1))
      for (a in LV) cells <- c(cells, sprintf("$%.3f$", mean(pv <= a, na.rm = TRUE)))
    }
    lines <- c(lines, sprintf("$%d$ & $%d$ & %s \\\\", T, N, paste(cells, collapse = " & ")))
  }
}
writeLines(lines, file.path(OUT, "tab_ndfcal.tex"))

# ---- Section 5 tables ---------------------------------------------
app <- read.csv(file.path(SIM, "results", "application", "application.csv"), stringsAsFactors = FALSE)
h <- app[app$r0 == 0.15 & app$criterion == "BIC" & app$normalization == "temp on sl", ]
sym <- c(FIEG = "\\textit{FIEG}(r_0)", BIEG = "\\textit{BIEG}(r_0,1)", GIEG = "\\textit{GIEG}(r_0)")
lines <- c()
for (case in c("c", "ct")) {
  panel <- if (case == "c") "A" else "B"
  lines <- c(lines,
    "\\midrule",
    sprintf("\\multicolumn{7}{l}{\\textit{Panel %s: Case $\\mathrm{%s}$}} \\\\", panel, case),
    "\\midrule")
  for (tst in TESTS) {
    r <- h[h$case == case & h$test == tst, ][1, ]
    lines <- c(lines, sprintf("$%s$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & $%.3f$ & %d--%d \\\\",
               sym[[tst]], r$stat, r$pvalue, r$cv1, r$cv5, r$cv10, r$win_start, r$win_end))
  }
}
writeLines(lines, file.path(OUT, "tab_app.tex"))

# robustness, case ct
lines <- c()
for (r0 in c(0.10, 0.15, 0.20)) for (crit in c("BIC", "AIC")) for (nrm in c("temp on sl", "sl on temp")) {
  cells <- c()
  for (tst in TESTS) {
    r <- app[app$case == "ct" & abs(app$r0 - r0) < 1e-9 & app$criterion == crit & app$normalization == nrm & app$test == tst, ][1, ]
    cells <- c(cells, sprintf("$%.3f$ & $%.3f$ & %d--%d", r$stat, r$pvalue, r$win_start, r$win_end))
  }
  ylab <- if (nrm == "temp on sl") "$y=$ temp." else "$y=$ s.\\,l."
  lines <- c(lines, sprintf("$%.2f$ & %s & %s & %s \\\\", r0, crit, ylab, paste(cells, collapse = " & ")))
}
writeLines(lines, file.path(OUT, "tab_robust.tex"))
cat("wrote tables to", OUT, "\n")

# ---- console: numbers for the prose -------------------------------
cat("\nNDF grid:", length(unique(ndf$prob)), "points\n")
f140 <- file.path(SIM, "results", "stage2", "stage2_case-c_T-140.csv")
f96  <- file.path(SIM, "results", "stage2", "stage2_case-c_T-96.csv")
if (file.exists(f140)) { cat("stage2 T=140:\n"); print(read.csv(f140)) }
if (file.exists(f96))  { cat("\nstage2 T=96:\n");  print(read.csv(f96)) }
cat("\n\nMCSE ranges, N=1:\n")
for (case in c("c", "ct")) for (tst in TESTS) {
  s <- q[q$case == case & q$N == 1 & q$test == tst & round(q$prob, 5) %in% LV, ]
  if (nrow(s) == 0) next
  cat(sprintf("  case %-2s %-4s mcse median %.3f max %.3f\n", case, tst,
              stats::median(s$mcse, na.rm = TRUE), max(s$mcse, na.rm = TRUE)))
}
