test_that("test_ss_fixedtrace() soundly doesn't reject for simulation of single sample from null", {
  set.seed(13131)
  Y <- rsymm_norm(200, diag(c(3,2,1)/6)) #50 matrices is too small
  Y <- lapply(Y, function(m) {m/sum(diag(m))})
  res <- test_ss_fixedtrace(Y, c(3,2,1)/6, 100)
  expect_gt(res$pval, 0.2)
  #hopefully even more accurate with correct evectors supplied, but doesn't seem to be
  res2 <- test_ss_fixedtrace(Y, c(3,2,1)/6, 100, evecs = diag(1, nrow = 3))
  expect_gt(res2$pval, 0.2)
})

test_that("test_ss_fixedtrace() reject for single sample with wrong eval", {
  set.seed(136)
  Y <- rsymm_norm(1000, diag(c(3,2,1)/6))
  Y <- lapply(Y, function(m) {m/sum(diag(m))})
  expect_warning(res <- test_ss_fixedtrace(Y, c(1,-1,1)/10, 100))
  expect_equal(res$pval, 0)
  
  badevals <- c(1,1,1)
  badevals <- badevals/sum(badevals)
  res <- test_ss_fixedtrace(Y, badevals, 100)
  expect_lt(res$pval, 0.05)
})

test_that("a simple multisample null situation ???", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) {m/sum(diag(m))})
    Y
  }, simplify = FALSE)
  
  stat_ss_fixedtrace(Ysamples[[1]], evals = c(3,2,1)/6)
  stat_ss_fixedtrace(Ysamples[[1]], evals = c(3,2,1)/6, evecs = diag(1, nrow = 3))
})


test_that("hasfixedtrace() gives TRUE or FALSE values", {
  set.seed(13)
  const <- 1.45
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) {const * m/sum(diag(m))})
    Y
    }, simplify = FALSE)

  expect_true(hasfixedtrace(Ysamples))
  expect_true(hasfixedtrace(Ysamples[[1]]))
  
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  
  expect_false(hasfixedtrace(Ysamples))
  expect_false(hasfixedtrace(Ysamples[[1]]))
})
