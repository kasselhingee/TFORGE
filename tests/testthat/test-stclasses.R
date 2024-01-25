test_that("as.mstorsst()",{
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  expect_silent(as.mstorsst(Ysamples))
  expect_silent(as.mstorsst(Ysamples[1]))
  
  expect_error(as.sst(Ysamples))
  expect_error(as.sst(matrix(1, 2, 5)), "number")
})

test_that("new generics work for sst", {
  set.seed(134)
  Ysample <- rsymm_norm(50, diag(c(3,2,1)))
  Ysample_sst <- as.sst(Ysample)
  expect_equal(Ysample_sst[1], Ysample[1])
  expect_equal(Ysample_sst[1:3], Ysample[1:3])
  boolindex <- runif(50) < 0.3
  expect_equal(Ysample_sst[boolindex], Ysample[boolindex])
  
  expect_equal(Ysample_sst[[1]], Ysample[[1]])
  expect_equal(Ysample_sst[[1:2]], Ysample[[1:2]])
})
