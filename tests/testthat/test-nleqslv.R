test_that("nleqslv error conditions work", {
  xstart <- c(1, 2, 3)
  z <- nleqslv(xstart, common_test_f, common_test_jac, method="Newton",
               control=list(trace=0, allowSingular=TRUE))
  expect_equal(z$x, common_test_xsol, tolerance = 1E-6)

  expect_error(z <- nleqslv(xstart, common_test_f, "Jac", method="Newton"))

  expect_error(z <- nleqslv(xstart, common_test_f, common_test_jac, method="Newton",
                            control=list(delta = as.Date("2026-01-01"))))
})
