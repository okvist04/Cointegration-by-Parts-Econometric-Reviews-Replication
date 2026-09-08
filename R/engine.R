# =====================================================================
# engine.R
#
# R-level wrapper around cbp_engine.cpp (the Rcpp port of fast_engine.jl).
# Compiles the C++ engine once via Rcpp::sourceCpp() and exposes an R API
# matching the Julia names/semantics as closely as R's conventions allow.
#
# Deterministic cases follow Equation (3.2) of the paper:
#   "n"   no deterministic term    (case_code 0)
#   "c"   constant                 (case_code 1)
#   "ct"  constant and linear trend (case_code 2)
#
# Correctness of the C++ engine itself is checked in verify_engine.R
# against a naive brute-force reference -- run that before trusting any
# number produced here.
# =====================================================================

library(Rcpp)

# Every script in this repository is run with the repository root as the
# working directory (see README.md), exactly as the Julia originals used
# `include(joinpath(@__DIR__, "fast_engine.jl"))` relative to their own
# location -- so this path is simply relative to the repo root.
sourceCpp("R/cbp_engine.cpp")

case_to_code <- function(case) {
  switch(case, n = 0L, c = 1L, ct = 2L, stop("case must be 'n', 'c', or 'ct'"))
}

det_dim <- function(case) {
  switch(case, n = 0L, c = 1L, ct = 2L, stop("case must be 'n', 'c', or 'ct'"))
}

pmax_rule <- function(n, mult = 12.0) pmax_rule_cpp(as.integer(n), mult)

min_window <- function(T, r0, floor_obs = 0L) min_window_cpp(as.integer(T), r0, as.integer(floor_obs))

# ---------------------------------------------------------------------
# Data generating processes (thin wrappers; RNG is R's own via
# set.seed(), so ordinary R reproducibility applies)
# ---------------------------------------------------------------------

dgp_null <- function(T, N, drift = 0.0) dgp_null_cpp(as.integer(T), as.integer(N), drift)

# ---------------------------------------------------------------------
# Parallel replication helper
#
# Julia's parallel_reps() hands replications to a pool of long-lived
# threads sharing a preallocated workspace. R's parallel::mclapply()
# instead forks one process per scheduled chunk; each replication here
# allocates its own (small, per-call) C++ workspace via rep_*_cpp(), so
# there is no persistent workspace to share, and fork overhead is
# negligible next to an O(T^2) scan for realistic T.
#
# Reproducibility: mclapply() with mc.set.seed = TRUE (the default)
# gives each forked child its own well-separated RNG stream, so results
# are reproducible given the same set.seed() call before the parallel
# section, but are NOT bitwise identical to a different mc.cores count
# (a limitation of fork-based parallel RNG in general, not specific to
# this code). For a fully deterministic single-machine record, run with
# mc.cores = 1.
# ---------------------------------------------------------------------

parallel_reps <- function(f, nrep, mc.cores = getOption("mc.cores", max(1L, parallel::detectCores() - 1L))) {
  if (.Platform$OS.type == "windows" || mc.cores <= 1) {
    lapply(seq_len(nrep), f)
  } else {
    parallel::mclapply(seq_len(nrep), f, mc.cores = mc.cores, mc.set.seed = TRUE)
  }
}

# ---------------------------------------------------------------------
# Collect a list-of-lists (one per replication) into a data.table-style
# named list of vectors -- the R equivalent of Julia's preallocated
# per-statistic vectors filled in the replication loop.
# ---------------------------------------------------------------------

collect_reps <- function(reps_list, fields) {
  out <- lapply(fields, function(fl) {
    vapply(reps_list, function(r) if (is.null(r) || is.null(r[[fl]])) NA_real_ else as.numeric(r[[fl]]), numeric(1))
  })
  names(out) <- fields
  out
}
