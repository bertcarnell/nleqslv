# Copyright (c) Rob Carnell 2026

test_that("Chebyquad solutions converge for n = 1:7,9 using default method", {
  chebyquad <- function(x) {
    n <- length(x)
    y <- numeric(n)

    for (j in 1:n) {
      t1 <- 1.0
      t2 <- 2.0 * x[j] - 1.0
      tmp <- 2.0 * t2

      for (i in 1:n) {
        y[i] <- y[i] + t2
        t3 <- tmp * t2 - t1
        t1 <- t2
        t2 <- t3
      }
    }

    y <- y / n

    for (i in 1:n) {
      if (i %% 2 == 0) {
        y[i] <- y[i] + 1.0 / (i * i - 1)
      }
    }

    y
  }

  chebyinit <- function(n) (1:n) / (n + 1)

  for (n in c(1:7, 9)) {
    xstart <- chebyinit(n)

    zz <- nleqslv(
      xstart, chebyquad,
      global = "dbldog",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        trace = 0,
        btol = 0.01,
        delta = -2
      )
    )

    expect_true(all(abs(zz$fvec) <= 1e-8))
  }
})


test_that("Chebyquad solutions converge for n = 1:7,9 using Newton method", {
  skip_if_not_installed("nleqslv")

  chebyquad <- function(x) {
    n <- length(x)
    y <- numeric(n)

    for (j in 1:n) {
      t1 <- 1.0
      t2 <- 2.0 * x[j] - 1.0
      tmp <- 2.0 * t2

      for (i in 1:n) {
        y[i] <- y[i] + t2
        t3 <- tmp * t2 - t1
        t1 <- t2
        t2 <- t3
      }
    }

    y <- y / n

    for (i in 1:n) {
      if (i %% 2 == 0) {
        y[i] <- y[i] + 1.0 / (i * i - 1)
      }
    }

    y
  }

  chebyinit <- function(n) (1:n) / (n + 1)

  for (n in c(1:7, 9)) {
    xstart <- chebyinit(n)

    zz <- nleqslv(
      xstart, chebyquad,
      global = "dbldog",
      method = "Newton",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        trace = 0,
        btol = 0.01,
        delta = -2
      )
    )

    expect_true(all(abs(zz$fvec) <= 1e-8))
  }
})
