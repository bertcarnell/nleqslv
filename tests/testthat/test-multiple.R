# Copyright (c) Rob Carnell 2026

test_that("multiple functions with nleqslv", {
  func_list <- list(brown, brown, dcbval, dciequ, helval, pwlsng, rosbrk,
                    vardim, watson, watson, wood, pwlbsc, dslnex,
                    f_nocedal)
  start_list <- list(brown_xstart(10), brown_xstart(20), dcbval_xstart(10),
                     dciequ_xstart(10), helval_xstart, pwlsng_xstart,
                     rosbrk_xstart, vardiminit(10), watson_xstart(6),
                     watson_xstart(9), wood_xstart, pwlbsc_xstart,
                     dslnex_xstart, nocedal_xstart)

  expect_equal(length(func_list), length(start_list))

  for (i in seq_along(func_list)) {
    soln <- nleqslv(start_list[[i]], func_list[[i]])
    expect_equal(soln$termcd, 1)
    expect_equal(soln$message, expectedMessage1)
    expect_true(all(abs(soln$fvec) <= 1e-7))
  }
})

test_that("multiple functions with testnslv", {
  func_list <- list(brown, dcbval, dciequ, helval, pwlsng, rosbrk,
                    vardim, watson, watson, wood, pwlbsc)
  start_list <- list(brown_xstart(20), dcbval_xstart(10), dciequ_xstart(10),
                     helval_xstart, pwlsng_xstart, rosbrk_xstart, vardiminit(10),
                     watson_xstart(6), watson_xstart(9), wood_xstart, pwlbsc_xstart)

  expect_equal(length(func_list), length(start_list))

  for (i in seq_along(func_list)) {
    soln <- testnslv(start_list[[i]], func_list[[i]])
    expect_true(inherits(soln, "test.nleqslv"))
    expect_true(is.data.frame(soln$out))
    expect_true(all(soln$out$termcd %in% c(1,4,5,6)))
    expect_true(all(soln$out$Iter < 92))
  }
})

test_that("vardim", {
  soln <- nleqslv(vardiminit(10), vardim)
  expect_equal(vardimsol(10), soln$x)

  soln <- nleqslv(vardiminit(5), vardim)
  expect_equal(vardimsol(5), soln$x)
})

test_that("with jacobian", {
  func_list <- list(pwlbsc)
  start_list <- list(pwlbsc_xstart)
  jac_list <- list(pwlbscjac)

  expect_equal(length(func_list), length(start_list))

  for (i in seq_along(func_list)) {
    soln <- nleqslv(start_list[[i]], func_list[[i]], jac_list[[i]])
    expect_equal(soln$termcd, 1)
    expect_equal(soln$message, expectedMessage1)
    expect_true(all(abs(soln$fvec) <= 1e-7))
  }
})

test_that("dslnex solutions", {
  z <- nleqslv(dslnex_xstart, dslnex, jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, global="none", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE, global="none")
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  dslnex_xstart <- c(2, .5)

  z <- nleqslv(dslnex_xstart, dslnex, jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  z <- nleqslv(dslnex_xstart, dslnex, method="Newton", jacobian=TRUE)
  expect_equal(jacdsln(z$x), z$jac, tolerance = 0.1, ignore_attr = TRUE)

  set.seed(11)

  dslnex_xstart <- matrix(runif(50, -3, 3), ncol=2)
  temp <- searchZeros(dslnex_xstart, dslnex)
  # two solutions of two values
  expect_equal(nrow(temp$x), 2)
  expect_equal(ncol(temp$x), 2)
  expect_true(all(apply(temp$x, 1, dslnex) < 1E-9))

})

test_that("nocedal", {
  z <- nleqslv(nocedal_xstart, f_nocedal, global="none")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  z <- nleqslv(nocedal_xstart, f_nocedal, method="Newton")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  z <- nleqslv(nocedal_xstart, f_nocedal, jac_nocedal, method="Newton")
  expect_equal(z$fvec, c(0, 0), tolerance = 1E-9)
  expect_equal(jac_nocedal(c(0,1)), jac_nocedal(z$x), tolerance = 1E-9)
})
