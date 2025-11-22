test_that("stat_ss1() on single sample from NULL is consistent with chisq", {
  set.seed(130)
  vals <- replicate(100, {
    Y <- rsymm_norm(30, diag(c(3,2,1)))
    Y <- normL2evals_sst(Y)
    stat_ss1(Y, evals = c(3,2,1)/sqrt(sum(c(3,2,1)^2)))
    })
  
  # qqplot(vals, y = rchisq(1E6, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1() on multiple NULL samples is consistent with chisq", {
  set.seed(13)
  vals <- replicate(100, {
  Ysamples <- lapply(c(2000,100,100,100), function(n){
    Y <- rsymm_norm(n, diag(c(3,2,1)))
    Y <- normL2evals_sst(Y)
    Y
  })
  stat_ss1(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (4-1)*2))
  res <- ks.test(vals, "pchisq", df = (4-1)*2)
  expect_gt(res$p.value, 0.2)
})

test_that("test_ss1() uniform pval on NULL TFORGE_fsm", {
  set.seed(1333)
  pvals <- replicate(ifelse(fast_check_on(), 20, 100), {
    Y <- rsymm_norm(30, diag(c(3,2,1)), sigma = diag(1, 6))
    Y <- normL2evals_sst(Y)
    res <- test_ss1(Y, normalise_ss1_vec(c(3,2,1)), ifelse(fast_check_on(), 20, 100), maxit = ifelse(fast_check_on(), 10, 1000))
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("chisq: test_ss1() uniform pval on NULL TFORGE_fsm", {
  set.seed(1333)
  pvals <- replicate(100, {
    Y <- rsymm_norm(30, diag(c(3,2,1)), sigma = diag(1, 6))
    Y <- normL2evals_sst(Y)
    res <- test_ss1(Y, normalise_ss1_vec(c(3,2,1)), "chisq", maxit = 1000)
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test_ss1() uniform pval on NULL TFORGE_kfsm", {
  set.seed(1333)
  pvals <- replicate(ifelse(fast_check_on(), 20, 100), {
    Ysamples <- replicate(2, {
      Y <- rsymm_norm(30, diag(c(3,2,1)))
      Y <- normL2evals_sst(Y)
    }, simplify = FALSE)
    res <- test_ss1(Ysamples, B = ifelse(fast_check_on(), 20, 100), maxit = ifelse(fast_check_on(), 10, 1000))
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("chisq: test_ss1() uniform pval on NULL TFORGE_kfsm", {
  set.seed(1333)
  pvals <- replicate(100, {
    Ysamples <- replicate(2, {
      Y <- rsymm_norm(30, diag(c(3,2,1)))
      Y <- normL2evals_sst(Y)
    }, simplify = FALSE)
    res <- test_ss1(Ysamples, B = "chisq", maxit = 1000)
    res$pval
  })
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test_ss1() reject for single sample with wrong eval", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(3,2,1)), sigma = diag(0.7, 6))
  Y <- normL2evals_sst(Y)
  res <- test_ss1(Y, normalise_ss1_vec(c(1,1,1)), ifelse(fast_check_on(), 10, 100), maxit = ifelse(fast_check_on(), 20, 100))
  expect_lt(res$pval, 0.05)
})

test_that("test_ss1() reject with warning for single sample with very wrong eval", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(3,2,1)), sigma = diag(0.7, 6))
  Y <- normL2evals_sst(Y)
  expect_warning(res <- test_ss1(Y, normalise_ss1_vec(c(-3,-2,-1)), 100, maxit = 100),
                 "convex hull")
  expect_equal(res$pval, 0)
})

test_that("test_ss1() reject for an TFORGE_kfsm", {
  set.seed(1333)
  Ysamples <- list(
    rsymm_norm(30, diag(c(3,2,1))),
    rsymm_norm(30, diag(c(6,2,1)))
  )
  Ysamples <- lapply(Ysamples, function(Y){
    normL2evals_sst(Y) #replace eigenvalues with normalised ones
  })
  res <- test_ss1(Ysamples, B = ifelse(fast_check_on(), 10, 100), maxit = ifelse(fast_check_on(), 10, 1000))
  expect_lt(res$pval, 0.05)
})

test_that("amaral2007Lemma1() produces correct result for a unit vector", {
  m <- runif(5, -1, 1)
  m <- m/sqrt(sum(m^2))
  
  A <- amaral2007Lemma1(m)
  expect_equal(drop(A %*% m), rep(0, 4))
  expect_equal(A %*% Conj(t(A)), diag(1, 4))
})

test_that("has_ss1 works for p=3", {
  set.seed(1333)
  Ysamples <- list(
    rsymm_norm(30, diag(c(3,2,1))),
    rsymm_norm(30, diag(c(6,2,1)))
  )
  expect_false(has_ss1(Ysamples))
  expect_error(test_ss1(Ysamples, B = 100))
  
  Ysamples <- lapply(Ysamples, function(Y){
    normL2evals_sst(Y) #replace eigenvalues with normalised ones
  })
  
  expect_true(has_ss1(Ysamples))
  
})
test_that("has_ss1 works for p=4", {
  set.seed(1333)
  Ysamples <- list(
    rsymm_norm(30, diag(c(4,3,2,1))),
    rsymm_norm(30, diag(c(7,6,2,1)))
  )
  expect_false(has_ss1(Ysamples))
  expect_error(test_ss1(Ysamples, B = 100))
  
  Ysamples <- lapply(Ysamples, function(Y){
    normL2evals_sst(Y) #replace eigenvalues with normalised ones
  })
  
  expect_true(has_ss1(Ysamples))
})
