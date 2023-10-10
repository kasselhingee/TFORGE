test_that("stat single sample has correct NULL distribution for projected trace", {
  set.seed(6514)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1) - 2))
    Y <- lapply(Y, projtrace) #this method of getting the correct trace seems to create narrower distributions than the normalising method
    stat_fixedtrace(Y, c(3,2,1) - 2)
  })

  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat single sample has WRONG NULL distribution for normtrace", {
  set.seed(6514)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1)/6), sigma = 0.05 * diag(1, 3*2))
    Y <- lapply(Y, normtrace) 
    stat_fixedtrace(Y, c(3,2,1)/6)
  })
  
  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_lt(res$p.value, 0.01)
})

test_that("stat on multi sample has correct NULL distribution", {
  set.seed(65141)
  vals <- replicate(100, {
    Ysamples <- replicate(5, {
      Y <- rsymm_norm(50, diag(c(3,2,1)))
      Y <- lapply(Y, projtrace)
      Y
    }, simplify = FALSE)
    stat_fixedtrace(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*2))
  res <- ks.test(vals, "pchisq", df = (5-1)*2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat on normed multi sample has correct NULL distribution", {
  set.seed(65141)
  vals <- replicate(100, {
    Ysamples <- replicate(5, {
      Y <- rsymm_norm(50, diag(c(3,2,1)))
      Y <- lapply(Y, normtrace)
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
    Y <- lapply(Y, projtrace)
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
      Y <- lapply(Y, projtrace)
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
  Y <- lapply(Y, projtrace)
  
  badevals <- c(1,1,-2)
  expect_error(res <- test_fixedtrace(Y, evals = badevals+1, B = 100))
  
  expect_warning(res <- test_fixedtrace(Y, badevals, B = 100, maxit = 100))
  expect_lt(res$pval, 0.05)
})

test_that("a multisample strongly non-null situation rejects", {
  set.seed(13)
  symm <- function(n, mn){
    Y <- rsymm_norm(n, mn)
    lapply(Y, projtrace)
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

test_that("projtrace() returns matrices that satisty hasfixedtrace", {
  set.seed(6514) 
  Y <- rsymm_norm(3, diag(c(3,2,-3)/6))
  Y <- lapply(Y, projtrace) #this method of getting the correct trace seems to create narrower distributions than the normalising method
  expect_true(hasfixedtrace(Y))
})


