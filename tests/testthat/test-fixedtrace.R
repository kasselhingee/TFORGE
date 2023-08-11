test_that("test_ss_fixedtrace() soundly doesn't reject for simulation of single sample from null", {
  set.seed(13131)
  Y <- rsymm_norm(50, diag(c(3,2,1)/6)) #50 matrices is too small
  Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)}) #this method of getting the correct trace seems to create narrower distributions than the normalising method
  stat_ss_fixedtrace(Y, c(3,2,1)/6)
  stat_ms_fixedtrace(Y, c(3,2,1)/6)
  stat_ss_fixedtrace(Y, c(3,2,1)/6, evecs = diag(1, nrow = 3))
  stat_ms_fixedtrace(Y, c(3,2,1)/6, evecs = diag(1, nrow = 3))
  res <- test_ss_fixedtrace(Y, c(3,2,1)/6, 100)
  expect_gt(res$pval, 0.2)
  #hopefully even more accurate with correct evectors supplied, but doesn't seem to be
  res2 <- test_ss_fixedtrace(Y, c(3,2,1)/6, 100, evecs = diag(1, nrow = 3))
  expect_gt(res2$pval, 0.2)
})

test_that("test_ss_fixedtrace() reject for single sample with wrong eval", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(3,2,1)/6), sigma = diag(rep(0.1, 6)))
  Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
  expect_warning(res <- test_ss_fixedtrace(Y, c(1,-1,1)/10, 100))
  expect_equal(res$pval, 0)
  
  badevals <- c(1,1,1)
  badevals <- badevals/sum(badevals)
  res <- test_ss_fixedtrace(Y, badevals, 100)
  expect_lt(res$pval, 0.05)
 
  #try with eigenvectors supplied
  res <- test_ss_fixedtrace(Y, badevals, 100, evecs = diag(1, 3))
  expect_lt(res$pval, 0.05)
})

test_that("a simple multisample null situation doesn't reject", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
    Y
  }, simplify = FALSE)
  
  stat_ms_fixedtrace(Ysamples)
  res <- test_ms_fixedtrace(Ysamples, 100)
  expect_gt(res$pval, 0.2)
})

test_that("a multisample strongly non-null situation rejects", {
  set.seed(13)
  symm <- function(n, mn){
    Y <- rsymm_norm(n, mn)
    lapply(Y, function(m) {diag(m) <- diag(m) - mean(diag(m)) + 1/3; return(m)})
  }
  Y1 <- symm(50, diag(c(3,2,1)))
  # lapply(Y1, function(m) sum(diag(m)))
  Ys <- list(Y1,
       symm(50, diag(c(1,1,1))))
  
  res <- test_ms_fixedtrace(Ys, 100)
  expect_lt(res$pval, 0.05)
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
