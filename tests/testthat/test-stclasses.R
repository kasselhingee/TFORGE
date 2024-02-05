test_that("as.mstorsst() on list of sst",{
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  expect_silent(as.mstorsst(Ysamples))
  expect_silent(as.mstorsst(Ysamples[1]))
  
  expect_error(as.sst(Ysamples))
  expect_error(as.sst(matrix(1, 2, 5)), "Number")
})

test_that("as.mstorsst() on list of list of matrices", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- apply(Y, 1, 
               function(vec) {m <- invvech(vec)},
               simplify = FALSE)
    Y
    }, simplify = FALSE)
  
  expect_s3_class(as.mstorsst(Ysamples), "mst")
})

test_that("matrix generics work for sst", {
  set.seed(134)
  Ysample <- matrix(runif(6*5), nrow = 5)
  Ysample_sst <- as.sst(Ysample)
  expect_equal(Ysample_sst[1], Ysample[1])
  expect_equal(Ysample_sst[1:3], Ysample[1:3])
  expect_equal(Ysample_sst[1,3], Ysample[1,3])
  boolindex <- runif(50) < 0.3
  expect_equal(Ysample_sst[boolindex], Ysample[boolindex])
  
  expect_equal(Ysample_sst[[3]], Ysample[[3]])
  expect_error(Ysample_sst[[1:2]], "more than one")
})
