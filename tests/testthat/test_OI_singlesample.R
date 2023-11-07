test_that("covOI gives correct matrices", {
  expect_equal(covOI(3, 2, 0, vectorisor = "vecd"),
               4 * diag(6))
  expect_equal(covOI(3, 1, 1/4, vectorisor = "vecd"),
               blockdiag(1 + diag(3), diag(3)))
  
  vechcovOI <- covOI(3, 1, 1/4, vectorisor = "vech")
  expect_equal(diag(vechcovOI), rep(2, 6))
  expect_equal(which(vechcovOI == 1, arr.ind = TRUE, useNames = FALSE),
  matrix(c(4,1,
           6,1,
           1,4,
           6,4,
           1,6,
           4,6), ncol = 2, byrow = TRUE), ignore_attr = TRUE)
})