test_that("makeSymmetric gives symmetric matrix for high dimension", {
  set.seed(34)
  m <- rsymm(1, diag(1:7))[[1]]
  magain <- makeSymmetric(m)
  expect_equal(m, magain)
  
  ess <- eigen_desc(m)
  expect_equal(m, makeSymmetric(ess$vectors %*% diag(ess$values) %*% t(ess$vectors)))
  
  m <- invvec(rnorm(7*7), nrow = 7)
  expect_error(makeSymmetric(m))
  expect_true(isSymmetric(makeSymmetric(m, tolerance = 1E2)))
})
