test_that("mmean is zero for simulated data", {
  Ysample <- rsymm(1000, matrix(0, 5, 5))
  expect_equal(mmean(Ysample), matrix(0, nrow = 5, ncol = 5), 
               tolerance = 2 * sqrt(2^2 / 12) / sqrt(1000))
})

test_that("merr returns input when input centred", {
  Ysample <- rsymm(10, matrix(0, 5, 5))
  Ysample <- lapply(Ysample, function(y) y - mmean(Ysample))
  expect_equal(mmean(Ysample), matrix(0, nrow = 5, ncol = 5))
  expect_equal(merr(Ysample), Ysample)
})

test_that("mcovar gets close to true", {
  skip("slow")
  Ysample <- rsymm(100000, matrix(0, 3, 3), sigma = diag( 3 * (3+1)/2))
  expect_silent(c0 <- mcovar(merr(Ysample)))
  expect_equal(c0, diag( 3 * (3+1)/2), tolerance = 0.06)
})

