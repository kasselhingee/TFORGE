test_that("regularellipse() point satisfy on equation", {
  a <- 3
  b <- 5
  expect_equal(drop(regularellipse(0, a = a, b = b)), c(x = a, y = 0))
  expect_equal(drop(regularellipse(pi/2, a = a, b = b)), c(x = 0, y = b))


  locs <- regularellipse(seq(-2*pi, 2*pi, length.out = 100), a = a, b = b)
  expect_equal(locs[, "x"]^2/(a^2) + locs[, "y"]^2/(b^2), rep(1, 100))
})

test_that("conf_fixedtrace() contains about 95% of sample points", {
  set.seed(345)
  ms <- rsymm_norm(n = 100, mean = diag(c(4,2,1)))
  ms <- normtrace(ms)
})
