# Copyright (c) Rob Carnell 2026

test_that("Trigonometric system solves correctly with global='dbldog'", {
  trig <- function(x) {
    n <- length(x)
    y <- cos(x)
    s <- sum(y)
    y <- n - s + seq_len(n) * (1 - y) - sin(x)
    y
  }

  trigjac <- function(x) {
    n <- length(x)
    J <- matrix(0, n, n)
    for (p in 1:n) {
      J[, p] <- sin(x[p])
      J[p, p] <- (p + 1) * sin(x[p]) - cos(x[p])
    }
    J
  }

  n <- 10
  xstart <- rep(1, n) / n

  znlm <- nleqslv(
    xstart,
    trig,
    global = "dbldog",
    control = list(trace = 0)
  )

  # Convergence checks
  expect_equal(znlm$termcd, 1)
  expect_equal(znlm$message, "Function criterion near zero")

  # Residual check
  expect_true(all(abs(znlm$fvec) <= 1e-8))
})
