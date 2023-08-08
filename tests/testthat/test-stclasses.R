test_that("as.mstorsst()",{
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  expect_silent(as.mstorsst(Ysamples))
  expect_silent(as.mstorsst(Ysamples[1]))
  
  expect_error(as.sst(Ysamples))
  expect_error(as.sst(matrix(1, 2, 3)))
})