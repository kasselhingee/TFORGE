test_that("stat_fixedtrace() single sample has correct NULL distribution", {
  set.seed(6514)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1)/6))
    Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)}) #this method of getting the correct trace seems to create narrower distributions than the normalising method
    stat_fixedtrace(Y, c(3,2,1)/6)
  })

  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_fixedtrace() multi sample has correct NULL distribution", {
  set.seed(65141)
  vals <- replicate(100, {
    Ysamples <- replicate(5, {
      Y <- rsymm_norm(50, diag(c(3,2,1)))
      Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
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
    Y <- rsymm_norm(50, diag(c(3,2,1)/6), sigma = diag(rep(0.1, 6)))
    Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
    res <- test_fixedtrace(Y, c(3,2,1)/6, 100, maxit = 100)
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
      Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
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
  Y <- rsymm_norm(50, diag(c(3,2,1)/6), sigma = diag(rep(0.1, 6)))
  Y <- lapply(Y, function(m) {m[1,1] <- 1 - sum(diag(m)[-1]); return(m)})
  expect_warning(res <- test_fixedtrace(Y, evals = c(1,-1,1)/10, B = 100))
  expect_equal(res$pval, 0)
  
  badevals <- c(1,1,1)
  badevals <- badevals/sum(badevals)
  res <- test_fixedtrace(Y, badevals, B = 100, maxit = 100)
  expect_lt(res$pval, 0.05)
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
