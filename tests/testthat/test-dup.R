test_that("dup() gives matrix that works on vec() and vech() results", {
  m <- matrix(c(1,2,3,2,5,6,3,6,9), byrow = FALSE, nrow = 3)
  expect_equal(drop(dup(3) %*% vech(m)), vec(m))
})
