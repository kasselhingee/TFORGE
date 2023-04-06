test_that("vecd() works", {
  expect_equal(vecd(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1, 5, 9, sqrt(2) * c(2,3,6)))
})

test_that("S_mcovar gets close to true", {
  skip("slow")
  Ysample <- rsymm(100000, matrix(0, 3, 3), sigma = diag( 3 * (3+1)/2))
  expect_silent(c0 <- S_mcovar(merr(Ysample)))
  expect_equal(c0, diag(c(rep(1, 3), rep(2, 3))), tolerance = 0.06)
})

test_that("makeblockdiagonal works", {
  A <- matrix(1, 3, 5)
  B <- matrix(2, 4, 6)
  expect_equal(blockdiag(A, B)[4,1], 0)
  expect_equal(blockdiag(A, B)[4,6], 2)
})

test_that("S_anv() gives exact distribution for 'OI' covariances", {
  stop()
})