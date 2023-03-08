testthat("mmean is zero for simulated data", {
  Ysample <- rsymm(1000, 5)
  expect_equal(mmean(Ysample), matrix(0, nrow = 5, ncol = 5), 
               tolerance = 2 * sqrt(2^2 / 12) / sqrt(1000))
})

estthat("merr returns input when input centred", {
  Ysample <- rsymm(10, 5)
  Ysample <- lapply(Ysample, function(y) y - mmean(Ysample))
  expect_equal(mmean(Ysample), matrix(0, nrow = 5, ncol = 5))
  expect_equal(merr(Ysample), Ysample)
})

estthat("mcovar runs", {
  Ysample <- rsymm(10, 5)
  expect_silent(c0 <- mcovar(merr(Ysample)))
  expect_equal(dim(c0), c(5*(5+1)/2, 5*(5+1)/2))
})

