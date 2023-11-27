test_that("regularellipse() point satisfy on equation", {
  a <- 3
  b <- 5
  expect_equal(drop(regularellipse(0, a = a, b = b)), c(x = a, y = 0))
  expect_equal(drop(regularellipse(pi/2, a = a, b = b)), c(x = 0, y = b))


  locs <- regularellipse(seq(-2*pi, 2*pi, length.out = 100), a = a, b = b)
  expect_equal(locs[, "x"]^2/(a^2) + locs[, "y"]^2/(b^2), rep(1, 100))
})

test_that("conf_fixedtrace() warns when ordered boundary is intersected", {
  set.seed(345)
  ms <- rsymm_norm(n = 10, mean = diag(c(4,2,1)))
  ms <- normtrace(ms)
  expect_warning(cr <- conf_fixedtrace(ms, alpha = 0.05, B = 1000, npts = 1000),
                 "boundary")
})
test_that("conf_fixedtrace() boundary plots nicely", {
  set.seed(345)
  ms <- rsymm_norm(n = 50, mean = diag(c(0.5,0.3,0.2)), sigma = diag(0.02, 6))
  ms <- normtrace(ms)
  cr <- conf_fixedtrace(ms, alpha = 0.05, B = 1000, npts = 1000)
  points(cr$boundary[, 1], cr$boundary[, 2])
  abline(a = 0, b = 1)
  abline(a = 1/2, b = -1/2)
  samplemeanevals <- t(replicate(100, {
    ms <- rsymm_norm(n = 50, mean = diag(c(4,2,1)))
    ms <- normtrace(ms)
    evals <- eigen_desc(mmean(ms))$values
  }))
  points(cr$est[1], cr$est[2], pch = 2)
  points(samplemeanevals[, 1], samplemeanevals[, 2], pch = 3)
  plot(samplemeanevals[, 1], samplemeanevals[, 2], pch = 3)
})

test_that("conf_fixedtrace() contains population mean about 95% of the time", {
  set.seed(345)
  popmeanincr <- pbapply::pbreplicate(500, {
     ms <- rsymm_norm(n = 30, mean = diag(c(4,2,1)))
     ms <- normtrace(ms)
     cr <- suppressWarnings(conf_fixedtrace(ms, alpha = 0.05, B = 1000, npts = 1000))
     cr$inregion(c(4,2,1)/7)})
  expect_equal(mean(popmeanincr), 0.95, tolerance = 5E-2)
})

