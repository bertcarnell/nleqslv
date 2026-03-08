# Copyright (c) Rob Carnell 2026

# Broyden banded function
brdban <- function(x) {
  ml <- 5
  mu <- 1
  n <- length(x)
  y <- numeric(n)

  for (k in 1:n) {
    k1 <- max(1, k - ml)
    k2 <- min(n, k + mu)

    temp <- 0
    for (j in k1:k2) {
      if (j != k) {
        temp <- temp + x[j] * (1 + x[j])
      }
    }

    y[k] <- x[k] * (2 + 5 * x[k]^2) + 1 - temp
  }
  y
}

test_that("Known solution produces correct function values", {
  xsol <- c(
    -0.42830, -0.47660, -0.51965, -0.55810, -0.59251,
    -0.62450, -0.62324, -0.62139, -0.62045, -0.58647
  )

  fsol <- brdban(xsol)

  expect_true(all(abs(fsol) < 1e-4))
})

test_that("nleqslv converges from xstart = -1", {
  n <- 10
  xstart <- rep(-1, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog", method = "Newton",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_true(all(abs(znlq$fvec) <= 1e-7))
})

test_that("nleqslv converges from xstart = -2 with Newton method", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog", method = "Newton",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_true(all(abs(znlq$fvec) <= 1e-8))
})

test_that("nleqslv converges from xstart = -2 without specifying method", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_true(all(abs(znlq$fvec) <= 1e-7))
})

test_that("Repeat convergence test for xstart = -2", {
  n <- 10
  xstart <- rep(-2, n)

  znlq <- nleqslv(
    xstart, brdban,
    global = "dbldog",
    control = list(trace = 0, ftol = 1e-8, xtol = 1e-8, btol = 1e-2, delta = -1)
  )

  expect_equal(znlq$termcd, 1)
  expect_true(all(abs(znlq$fvec) <= 1e-8))
})
