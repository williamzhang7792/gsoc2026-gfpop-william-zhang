# PoissonFPOP.R
# quick sanity check: compile solver, compare output with Segmentor3IsBack

library(Rcpp)
library(Segmentor3IsBack)

sourceCpp("PoissonFPOP.cpp")

# pick the best K via penalized cost: cost[k] + penalty*(k-1)
select_segmentor_model <- function(seg, penalty) {
  costs <- seg@likelihood
  K_values <- seq_len(nrow(costs))
  penalized <- costs[,1] + penalty * (K_values - 1)
  best_K <- which.min(penalized)
  breaks_raw <- getBreaks(seg)[best_K, ]
  params_raw <- getParameters(seg)[best_K, ]
  # drop trailing zeros (unused slots)
  breaks <- breaks_raw[seq_len(best_K)]
  params <- params_raw[seq_len(best_K)]
  list(K = best_K, breaks = breaks, means = params,
       cost = costs[best_K, 1], penalized_cost = penalized[best_K])
}

# test data with clear changepoints
set.seed(42)
data_vec <- as.integer(c(rpois(30, 5), rpois(30, 20), rpois(30, 5)))

cat("Data (first 10):", head(data_vec, 10), "...\n")
cat("Data length:", length(data_vec), "\n\n")

# ground truth
seg <- Segmentor(data_vec, model = 1, Kmax = 30)

# compare for a range of penalties
penalties <- c(1, 5, 10, 50, 100, 500)

for (pen in penalties) {
  fpop_result <- PoissonFPOP(data_vec, pen)
  seg_result <- select_segmentor_model(seg, pen)

  cat(sprintf("penalty = %6.1f | FPOP segments: %d | Segmentor segments: %d",
              pen, fpop_result$n.segments, seg_result$K))

  # convert FPOP 0-based ends to Segmentor 1-based breaks
  K <- fpop_result$n.segments
  fpop_breaks <- c(fpop_result$end[-1] + 1L, length(data_vec))
  seg_breaks <- seg_result$breaks

  if (length(fpop_breaks) == length(seg_breaks) &&
      all(fpop_breaks == seg_breaks)) {
    cat(" | breaks: MATCH")
  } else {
    cat(" | breaks: MISMATCH!")
    cat(sprintf("\n  FPOP: %s", paste(fpop_breaks, collapse=" ")))
    cat(sprintf("\n  Seg:  %s", paste(seg_breaks, collapse=" ")))
  }

  # means
  fpop_means <- round(fpop_result$mean, 4)
  seg_means <- round(seg_result$means, 4)
  if (length(fpop_means) == length(seg_means) &&
      all(abs(fpop_means - seg_means) < 0.01)) {
    cat(" | means: MATCH\n")
  } else {
    cat(" | means: MISMATCH!\n")
    cat("  FPOP means:", fpop_means, "\n")
    cat("  Seg  means:", seg_means, "\n")
  }
}

cat("\nDone.\n")
