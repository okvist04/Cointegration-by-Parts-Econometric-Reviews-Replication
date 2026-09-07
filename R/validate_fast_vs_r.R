# ============================================================
# validate_fast_vs_r.R
#
# Run this BEFORE trusting the fast (Rcpp) path for the real job.
#
# It compares minstat_one_rep_null() (fast, C++-backed) against
# minstat_one_rep_null_r() (original pure-R logic) on identical
# simulated data, across a range of T, n_regressors, and scan_type,
# and reports the max absolute difference in the returned min-t
# statistic. These should match to within floating-point tolerance
# (~1e-8 or smaller) -- anything larger means something in the C++
# translation needs a closer look before you rely on it.
#
# Deliberately kept small (small T, few reps) since "all" in pure R
# is exactly the slow path we're trying to avoid for the real job.
# ============================================================

source("R/generate_cv.R")

OUTPUT_DIR <- "output/cv_output"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

set.seed(1)

T_values <- c(60, 120, 250)
n_regressors_values <- c(1, 3, 6)
scan_types <- c("fw", "bw", "all")
n_rep_check <- 15   # small on purpose -- this exercises the pure-R
                     # "all" path, which is the slow one

rows <- list()
idx <- 1

total_comparisons <- length(T_values) * length(n_regressors_values) *
  length(scan_types) * n_rep_check

pb <- txtProgressBar(min = 0, max = total_comparisons, style = 3)

for (T in T_values) {
  for (k in n_regressors_values) {
    for (st in scan_types) {

      for (r in seq_len(n_rep_check)) {

        s <- 1000 + idx  # shared seed -> identical simulated data

        fast_val <- minstat_one_rep_null(
          T = T, n_regressors = k, scan_type = st, seed = s
        )

        r_val <- minstat_one_rep_null_r(
          T = T, n_regressors = k, scan_type = st, seed = s
        )

        rows[[idx]] <- data.table(
          T = T, n_regressors = k, scan_type = st, rep = r,
          fast = fast_val, r = r_val,
          abs_diff = abs(fast_val - r_val)
        )

        setTxtProgressBar(pb, idx)

        idx <- idx + 1
      }
    }
  }
}

close(pb)

out <- rbindlist(rows)

cat("====================================\n")
cat("Validation: fast (C++) vs pure R\n")
cat("====================================\n")
cat("Comparisons run :", nrow(out), "\n")
cat("Max abs diff    :", max(out$abs_diff, na.rm = TRUE), "\n")
cat("Mean abs diff   :", mean(out$abs_diff, na.rm = TRUE), "\n")
cat("NA mismatches   :", sum(is.na(out$fast) != is.na(out$r)), "\n")

worst <- out[order(-abs_diff)][1:10]
cat("\nWorst 10 comparisons:\n")
print(worst)

fwrite(out, file.path(OUTPUT_DIR, "validation_fast_vs_r.csv"))

cat("\nFull comparison table written to", file.path(OUTPUT_DIR, "validation_fast_vs_r.csv"), "\n")
cat("If max abs diff is tiny (e.g. < 1e-6) and there are no NA\n")
cat("mismatches, the fast path is safe to use for the real run.\n")
