# test-NormalFPOP.R
# Validate NormalFPOP (regularized isotonic regression, Normal loss)
# against isoreg() (penalty=0) and fpop::Fpop() (non-decreasing data).

library(testthat)
library(Rcpp)

sourceCpp("NormalFPOP.cpp")

# Helper: expand segment means to per-point fitted values
# result$end[k] = last index of segment k-1 (0-based), with end[1] = -1 sentinel
expand_means <- function(result, n) {
  fitted <- numeric(n)
  for (k in seq_len(result$n.segments)) {
    seg_start <- result$end[k] + 1L
    if (k < result$n.segments) {
      seg_stop <- result$end[k + 1L]
    } else {
      seg_stop <- n - 1L
    }
    idx <- (seg_start:seg_stop) + 1L
    fitted[idx] <- result$mean[k]
  }
  fitted
}

# -----------------------------------------------------------------------
# Test 0: Smoke test — strictly decreasing data, penalty=0
# -----------------------------------------------------------------------

test_that("strictly decreasing data gives flat mean (PAVA pooling)", {
  y <- c(5, 4, 3, 2, 1)
  result <- NormalFPOP(y, 0)
  fitted <- expand_means(result, length(y))
  expect_equal(fitted, rep(mean(y), length(y)), tolerance = 1e-6)
})

# -----------------------------------------------------------------------
# Test 1: penalty=0 vs isoreg()
# -----------------------------------------------------------------------

test_that("penalty=0 matches isoreg on random data", {
  set.seed(42)
  y <- c(5, 3, 4, 1, 2, 8, 7, 9)
  result <- NormalFPOP(y, 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "random 8-element vector")
})

test_that("penalty=0 matches isoreg on already sorted data", {
  y <- c(1, 2, 3, 4, 5)
  result <- NormalFPOP(y, 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "already non-decreasing")
})

test_that("penalty=0 matches isoreg on constant data", {
  y <- rep(3.0, 8)
  result <- NormalFPOP(y, 0)
  expect_equal(result$n.segments, 1L)
  expect_equal(result$mean[1], 3.0, tolerance = 1e-6)
})

test_that("penalty=0 matches isoreg on strictly decreasing data", {
  y <- c(10, 8, 6, 4, 2)
  result <- NormalFPOP(y, 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "strictly decreasing")
})

test_that("penalty=0 matches isoreg on longer random data", {
  set.seed(2026)
  y <- rnorm(50, mean = cumsum(rnorm(50, 0, 0.3)))
  result <- NormalFPOP(y, 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "50-point random walk")
})

# -----------------------------------------------------------------------
# Test 2: Non-decreasing data vs fpop::Fpop()
# -----------------------------------------------------------------------

if (requireNamespace("fpop", quietly = TRUE)) {

  test_that("non-decreasing data matches fpop::Fpop at several penalties", {
    set.seed(123)
    y <- c(rnorm(30, 2, 0.5), rnorm(30, 5, 0.5), rnorm(30, 10, 0.5))
    n <- length(y)

    for (pen in c(1, 5, 10, 50)) {
      our   <- NormalFPOP(y, pen)
      fpop_result <- fpop::Fpop(y, pen)

      our_fitted  <- expand_means(our, n)

      fpop_breaks <- fpop_result$t.est
      fpop_means  <- numeric(n)
      seg_start <- 1L
      for (b in fpop_breaks) {
        fpop_means[seg_start:b] <- mean(y[seg_start:b])
        seg_start <- b + 1L
      }

      our_means_nondecr <- all(diff(our$mean) >= -1e-8)
      expect_true(our_means_nondecr,
                  info = paste("means non-decreasing at penalty =", pen))

      expect_equal(our_fitted, fpop_means, tolerance = 1e-4,
                   info = paste("fitted values at penalty =", pen))
    }
  })

} else {
  message("fpop package not installed, skipping Fpop comparison tests")
}

# -----------------------------------------------------------------------
# Test 3: Edge cases
# -----------------------------------------------------------------------

test_that("n=1 gives trivial result", {
  result <- NormalFPOP(c(42.0), 5)
  expect_equal(result$n.segments, 1L)
  expect_equal(result$mean[1], 42.0)
  expect_equal(result$cost, 0.0)
})

test_that("n=2 works correctly", {
  y <- c(1.0, 5.0)
  # low penalty: 2 segments, means = 1 and 5
  r_low <- NormalFPOP(y, 0)
  expect_equal(r_low$n.segments, 2L)
  expect_equal(as.numeric(r_low$mean), c(1, 5), tolerance = 1e-6)

  # high penalty: 1 segment, mean = 3
  r_high <- NormalFPOP(y, 1e6)
  expect_equal(r_high$n.segments, 1L)
  expect_equal(r_high$mean[1], mean(y), tolerance = 1e-6)
})

test_that("means are always non-decreasing", {
  set.seed(7)
  y <- rnorm(30, mean = 0, sd = 5)
  for (pen in c(0, 1, 5, 20, 100)) {
    result <- NormalFPOP(y, pen)
    if (result$n.segments > 1) {
      expect_true(all(diff(result$mean) >= -1e-8),
                  info = paste("non-decreasing at penalty =", pen))
    }
  }
})

test_that("high penalty gives 1 segment at overall mean", {
  set.seed(99)
  y <- rnorm(40, mean = 5)
  result <- NormalFPOP(y, 1e8)
  expect_equal(result$n.segments, 1L)
  expect_equal(result$mean[1], mean(y), tolerance = 1e-6)
})

test_that("already isotonic data is not altered", {
  y <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  for (pen in c(0, 1, 5)) {
    result <- NormalFPOP(y, pen)
    fitted <- expand_means(result, length(y))
    iso_fitted <- as.numeric(isoreg(y)$yf)
    if (pen == 0) {
      expect_equal(fitted, iso_fitted, tolerance = 1e-6,
                   info = "isotonic data, penalty=0")
    }
    expect_true(all(diff(result$mean) >= -1e-8),
                info = paste("isotonic check at penalty =", pen))
  }
})

# -----------------------------------------------------------------------
# Test 4: Random fuzz — penalty=0 must match isoreg on many inputs
# -----------------------------------------------------------------------

test_that("penalty=0 matches isoreg on 200 random trials (N=10)", {
  set.seed(20260222)
  n_fail <- 0L
  for (trial in seq_len(200)) {
    y <- switch(((trial - 1L) %% 3L) + 1L,
      rnorm(10), round(runif(10, 0, 10)), rnorm(10, sd = 1e3))
    if (length(unique(y)) < 2L) next
    result <- NormalFPOP(y, 0)
    fitted <- expand_means(result, length(y))
    iso <- as.numeric(isoreg(y)$yf)
    if (max(abs(fitted - iso)) > 1e-5) n_fail <- n_fail + 1L
  }
  expect_equal(n_fail, 0L)
})

cat("\nAll tests passed!\n")
