test_that("descending order restarts handled - obsolete", {
  expect_output({out <- withCallingHandlers(
    descendingordererror(c(3,1,2)),
    est_evals_not_descending = function(e) {
      print(e$message)
      invokeRestart("use_NA")}
  )}, "3 1 2")
  expect_equal(out, expected = rep(NA_real_, 3))
  
  expect_silent(out2 <- withCallingHandlers(
    descendingordererror(c(3,1,2)),
    est_evals_not_descending = function(e) {
      invokeRestart("use_NA")}
  ))
  expect_equal(out2, expected = rep(NA_real_, 3))
  
  out <- withCallingHandlers(
    descendingordererror(c(3,1,2)),
    est_evals_not_descending = function(e) {invokeRestart("ignore")})
  expect_equal(out, expected = c(3,1,2))
  
  out <- withCallingHandlers(
    descendingordererror(c(3,1,2)),
    est_evals_not_descending = function(e) {invokeRestart("sort")})
  expect_equal(out, expected = c(3,2,1))
  
  out <- withCallingHandlers(
    descendingordererror(c(3,1,2)),
    est_evals_not_descending = function(e) {invokeRestart("use_value", "hello")})
  expect_equal(out, expected = "hello")
  
  expect_error(descendingordererror(c(3,1,2)))
})

test_that("descending order error activates in resampling", {
  #based on finding a situation in simstudy31
  set.seed(224)
  allsim <- list(
    rsymm_norm(15, mean = diag(c(4,2,1))),
    rsymm_norm(15, mean = diag(c(4,2,1)))
  )
  allsim <- lapply(allsim, normalise_trace)
  expect_warning(res <- test_fixedtrace(allsim, B = 10), "1 bootstrap")
  expect_gt(sum(grepl("not in descending order", res$nullt_messages)), 0)
})

test_that("singularity error activates in resampling", {
  #tiny sample to try to get error in covariance with resampling
  set.seed(3)
  s1 <- rsymm_norm(5, mean = diag(c(4,2,1)))
  allsim <- list( s1, s1 )
  allsim <- lapply(allsim, normalise_trace)
  expect_warning(res <- test_fixedtrace(allsim, B = 10), "bootstrap resamples")
  expect_gt(sum(grepl("singular", res$nullt_messages)), 0)
})

test_that("stat single sample has correct NULL distribution for projected trace", {
  set.seed(6514)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1) - 2))
    Y <- projtrace_sst(Y) #this method of getting the correct trace seems to create narrower distributions than the normalising method
    stat_fixedtrace(Y, c(3,2,1) - 2)
  })

  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat single sample has WRONG NULL distribution for normalise_trace", {
  set.seed(6514)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1)/6), sigma = 0.05 * diag(1, 3*2))
    Y <- normalise_trace(Y) 
    stat_fixedtrace(Y, c(3,2,1)/6)
  })
  
  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_lt(res$p.value, 0.01)
})

test_that("stat on multi sample has correct NULL distribution", {
  set.seed(65142)
  vals <- replicate(100, {
    Ysamples <- lapply(c(2000, 100, 100, 100), function(n){
      Y <- rsymm_norm(n, diag(c(3,2,1)))
      Y <- projtrace_sst(Y)
      Y
    })
    stat_fixedtrace(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (4-1)*2))
  res <- ks.test(vals, "pchisq", df = (4-1)*2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat on normed multi sample has correct NULL distribution", {
  set.seed(65141)
  vals <- replicate(100, {
    Ysamples <- replicate(5, {
      Y <- rsymm_norm(50, diag(c(3,2,1)))
      Y <- normalise_trace(Y)
      Y
    }, simplify = FALSE)
    stat_fixedtrace(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*2))
  res <- ks.test(vals, "pchisq", df = (5-1)*2)
  expect_gt(res$p.value, 0.2)
})


test_that("test of NULL has uniform p values for sst", {
  set.seed(1333)
  pvals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1) - 2), sigma = diag(rep(0.1, 6)))
    Y <- projtrace_sst(Y)
    res <- test_fixedtrace(Y, c(3,2,1) - 2, 100, maxit = 100)
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  expect_gt(suppressWarnings(ks.test(pvals, "punif")$p.value), 0.2)
})

test_that("test of NULL has uniform p values for mst", {
  set.seed(1333)
  pvals <- replicate(100, {
    Ysamples <- replicate(5, {
      Y <- rsymm_norm(50, diag(c(3,2,1)))
      Y <- projtrace_sst(Y)
      Y
    }, simplify = FALSE)
    res <- test_fixedtrace(Ysamples, B = 100, maxit = 100)
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  expect_gt(suppressWarnings(ks.test(pvals, "punif")$p.value), 0.2)
})

test_that("test rejects for single sample with wrong eval", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(1,0,-1)), sigma = diag(rep(0.1, 6)))
  Y <- projtrace_sst(Y)
  
  badevals <- c(1,1,-2)
  expect_error(res <- test_fixedtrace(Y, evals = badevals+1, B = 100))
  
  expect_warning(res <- test_fixedtrace(Y, badevals, B = 100, maxit = 100))
  expect_lt(res$pval, 0.05)
})

test_that("a multisample strongly non-null situation rejects", {
  set.seed(13)
  symm <- function(n, mn){
    Y <- rsymm_norm(n, mn)
    projtrace_sst(Y)
  }
  Y1 <- symm(50, diag(c(3,2,1)))
  # lapply(Y1, function(m) sum(diag(m)))
  Ys <- list(Y1,
       symm(50, diag(c(1,1,1))))
  
  res <- test_fixedtrace(Ys, B = 100, maxit = 100)
  expect_lt(res$pval, 0.05)
})


test_that("hasfixedtrace() gives TRUE or FALSE values", {
  set.seed(13)
  const <- 1.45
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- apply(Y, 1, 
               function(vec) {m <- invvech(vec); const * m/sum(diag(m))},
               simplify = FALSE)
    Y
    }, simplify = FALSE)
  
  expect_equal(sum(diag(Ysamples[[1]][[1]])), const)

  expect_true(hasfixedtrace(as.mstorsst(Ysamples)))
  expect_true(hasfixedtrace(as_fsm(Ysamples[[1]])))
  
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  
  expect_false(hasfixedtrace(Ysamples))
  expect_false(hasfixedtrace(Ysamples[[1]]))
})

test_that("projtrace() returns matrices that satisty hasfixedtrace", {
  set.seed(6514) 
  Y <- rsymm_norm(3, diag(c(3,2,-3)/6))
  Y <- projtrace_sst(Y) #this method of getting the correct trace seems to create narrower distributions than the normalising method
  expect_true(hasfixedtrace(Y))
})


