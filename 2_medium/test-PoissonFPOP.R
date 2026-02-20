# test-PoissonFPOP.R
# check that PoissonFPOP matches Segmentor3IsBack for various inputs/penalties

library(testthat)
library(Rcpp)
library(Segmentor3IsBack)

sourceCpp("PoissonFPOP.cpp")

# pick the best K for a given penalty via penalized cost
select_segmentor_model <- function(seg, penalty) {
  costs <- seg@likelihood
  K_values <- seq_len(nrow(costs))
  penalized <- costs[,1] + penalty * (K_values - 1)
  best_K <- which.min(penalized)
  breaks <- getBreaks(seg)[best_K, seq_len(best_K)]
  means  <- getParameters(seg)[best_K, seq_len(best_K)]
  list(K = best_K, breaks = breaks, means = means,
       cost = costs[best_K, 1])
}

# FPOP ends -> Segmentor-style 1-based breaks
fpop_to_breaks <- function(fpop_result, n) {
  c(fpop_result$end[-1] + 1L, n)
}

# --- 3-segment Poisson data ---

test_that("3-segment Poisson data matches Segmentor", {
  set.seed(42)
  data_vec <- as.integer(c(rpois(30, 5), rpois(30, 20), rpois(30, 5)))
  n <- length(data_vec)
  seg <- Segmentor(data_vec, model = 1, Kmax = 30)

  for (pen in c(5, 10, 50, 100)) {
    fpop <- PoissonFPOP(data_vec, pen)
    segm <- select_segmentor_model(seg, pen)
    fpop_breaks <- fpop_to_breaks(fpop, n)

    expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                     info = paste("n.segments at penalty =", pen))
    expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                     info = paste("breaks at penalty =", pen))
    expect_equal(as.numeric(fpop$mean), as.numeric(segm$means), tolerance = 1e-6,
                 info = paste("means at penalty =", pen))
  }
})

# --- single changepoint ---

test_that("single changepoint data matches Segmentor", {
  set.seed(123)
  data_vec <- as.integer(c(rpois(50, 3), rpois(50, 15)))
  n <- length(data_vec)
  seg <- Segmentor(data_vec, model = 1, Kmax = 20)

  for (pen in c(10, 50, 200)) {
    fpop <- PoissonFPOP(data_vec, pen)
    segm <- select_segmentor_model(seg, pen)
    fpop_breaks <- fpop_to_breaks(fpop, n)

    expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                     info = paste("n.segments at penalty =", pen))
    expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                     info = paste("breaks at penalty =", pen))
    expect_equal(as.numeric(fpop$mean), as.numeric(segm$means), tolerance = 1e-6,
                 info = paste("means at penalty =", pen))
  }
})

# --- constant rate (expect 1 segment at high penalty) ---

test_that("constant-rate data gives 1 segment at high penalty", {
  set.seed(7)
  data_vec <- as.integer(rpois(80, 10))
  n <- length(data_vec)
  seg <- Segmentor(data_vec, model = 1, Kmax = 10)

  fpop <- PoissonFPOP(data_vec, 500)
  segm <- select_segmentor_model(seg, 500)
  fpop_breaks <- fpop_to_breaks(fpop, n)

  expect_identical(as.integer(fpop$n.segments), 1L)
  expect_identical(as.integer(fpop$n.segments), as.integer(segm$K))
  expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks))
  expect_equal(as.numeric(fpop$mean), as.numeric(segm$means), tolerance = 1e-6)
})

# --- many changepoints, low penalty ---

test_that("many-changepoint data matches Segmentor at low penalty", {
  set.seed(99)
  data_vec <- as.integer(c(
    rpois(20, 2), rpois(20, 15), rpois(20, 4),
    rpois(20, 25), rpois(20, 3)
  ))
  n <- length(data_vec)
  seg <- Segmentor(data_vec, model = 1, Kmax = 30)

  for (pen in c(1, 5, 20, 100)) {
    fpop <- PoissonFPOP(data_vec, pen)
    segm <- select_segmentor_model(seg, pen)
    fpop_breaks <- fpop_to_breaks(fpop, n)

    expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                     info = paste("n.segments at penalty =", pen))
    expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                     info = paste("breaks at penalty =", pen))
    expect_equal(as.numeric(fpop$mean), as.numeric(segm$means), tolerance = 1e-6,
                 info = paste("means at penalty =", pen))
  }
})

# --- zeros in data (Poisson rate ~ 0.5 generates lots of zeros) ---

test_that("data containing zeros is handled correctly", {
  set.seed(2026)
  data_vec <- as.integer(c(rpois(30, 0.5), rpois(30, 10), rpois(30, 0.5)))
  n <- length(data_vec)
  seg <- Segmentor(data_vec, model = 1, Kmax = 30)

  for (pen in c(1, 5, 50)) {
    fpop <- PoissonFPOP(data_vec, pen)
    segm <- select_segmentor_model(seg, pen)
    fpop_breaks <- fpop_to_breaks(fpop, n)

    expect_identical(as.integer(fpop$n.segments), as.integer(segm$K),
                     info = paste("n.segments at penalty =", pen))
    expect_identical(as.integer(fpop_breaks), as.integer(segm$breaks),
                     info = paste("breaks at penalty =", pen))
    expect_equal(as.numeric(fpop$mean), as.numeric(segm$means), tolerance = 1e-6,
                 info = paste("means at penalty =", pen))
  }
})

# --- penalty=0 edge case ---

test_that("penalty=0 gives every point its own segment", {
  data_vec <- as.integer(c(3, 10, 3, 10, 3))
  fpop <- PoissonFPOP(data_vec, 0.0)
  expect_identical(as.integer(fpop$n.segments), length(data_vec))
  expect_equal(as.numeric(fpop$mean), as.numeric(data_vec))
})

# --- n=2 (smallest valid input) ---

test_that("n=2 works for both 1-segment and 2-segment outcomes", {
  data_vec <- as.integer(c(1, 100))
  # low penalty -> 2 segments
  fpop2 <- PoissonFPOP(data_vec, 0.0)
  expect_identical(as.integer(fpop2$n.segments), 2L)
  expect_equal(as.numeric(fpop2$mean), c(1, 100))
  # huge penalty -> 1 segment
  fpop1 <- PoissonFPOP(data_vec, 1e10)
  expect_identical(as.integer(fpop1$n.segments), 1L)
  expect_equal(as.numeric(fpop1$mean), mean(data_vec), tolerance = 1e-6)
})

# --- bad input ---

test_that("rejects negative data and all-identical data", {
  expect_error(PoissonFPOP(as.integer(c(-1, 2, 3)), 5), "Negative")
  expect_error(PoissonFPOP(as.integer(c(10, 10, 10)), 5), "identical")
})

cat("\nAll tests passed!\n")
