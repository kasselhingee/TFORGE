test_that("as_flat() on list of TFORGE_fsm",{
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  expect_silent(as_flat(Ysamples))
  expect_silent(as_flat(Ysamples[1]))
  
  expect_error(as_fsm(Ysamples))
  expect_error(as_fsm(matrix(1, 2, 5)), "Number")
})

test_that("as_flat() on list of list of matrices", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- apply(Y, 1, 
               function(vec) {m <- inv_vech(vec)},
               simplify = FALSE)
    Y
    }, simplify = FALSE)
  
  expect_s3_class(as_flat(Ysamples), "TFORGE_kfsm")
  
  Ysamples <- list({
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- apply(Y, 1, 
               function(vec) {m <- inv_vech(vec)},
               simplify = FALSE)
    Y
    }, {
    Y <- rsymm_norm(50, diag(c(3, 3,2,1)))
    Y <- apply(Y, 1, 
               function(vec) {m <- inv_vech(vec)},
               simplify = FALSE)
    Y})
  
  expect_error(as_flat(Ysamples), "different")
})

test_that("as_flat() on list of TFORGE_fsm-like matrices", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    class(Y) <- "matrix"
    Y
    }, simplify = FALSE)
  
  expect_s3_class(as_flat(Ysamples), "TFORGE_kfsm")
  
  Ysamples <- list({
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    class(Y) <- "matrix"
    Y
    }, {
    Y <- rsymm_norm(50, diag(c(3, 3,2,1)))
    class(Y) <- "matrix"
    Y})
  
  expect_error(as_flat(Ysamples), "different")
})


test_that("as_flat() on a matrix-like TFORGE_fsm", {
  set.seed(13)
  Y <- rsymm_norm(50, diag(c(3,2,1)))
  class(Y) <- "matrix"
  Ysst <- as_flat(Y)
  expect_s3_class(Ysst, "TFORGE_fsm")
  expect_equal(nrow(Ysst), nrow(Y))
  expect_equal(ncol(Ysst), ncol(Y))
})

test_that("as_fsm() error when bad lists of matrices passed", {
  expect_error(as_fsm(list(diag(3), "test")))
  expect_error(as_fsm(list(rbind(diag(3), rep(1, 3)), diag(4))))
  expect_error(as_fsm(list(diag(3), diag(4))))
  expect_error(as_fsm(list(diag(3), cbind(diag(3), rep(1, 3)))))
})

test_that("matrix generics work for TFORGE_fsm", {
  set.seed(134)
  Ysample <- matrix(runif(6*5), nrow = 5)
  Ysample_sst <- as_fsm(Ysample)
  expect_equal(Ysample_sst[1], Ysample[1])
  expect_equal(Ysample_sst[1:3], Ysample[1:3])
  expect_equal(Ysample_sst[1,3], Ysample[1,3])
  boolindex <- runif(50) < 0.3
  expect_equal(Ysample_sst[boolindex], Ysample[boolindex])
  
  expect_equal(Ysample_sst[[3]], Ysample[[3]])
  expect_error(Ysample_sst[[1:2]], "more than one")
})

test_that("as_flat() works on a list of matrices", {
  x <- rsymm_Schwartzman(10, diag(c(1,2,4)), drop(rWishart(1, 6, diag(6))))
  res <- as_flat(x)
  expect_s3_class(res, "TFORGE_fsm")
  expect_equal(dim(res), c(10, 6))
})
