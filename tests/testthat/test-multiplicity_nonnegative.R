test_that("rsymm_lognorm normalised has correct eigenvalues", {
  skip_if_fast_check() #very slow to test
  # simulating from a narrow distribution to keep within bounds
  sdlog <- 0.1
  meanlogval <- log(1/3) - sdlog^2/2
  set.seed(3)
  Ysample <- rsymm_lognorm(1E3, meanlog = diag(meanlogval, 3), sigmalog = diag(sdlog^2, 6))
  
  #without projection/normalisation, the mean is slightly too high
  expect_true(all(eigen_desc(mmean(Ysample))$values > 1/3))
  
  # normalisation of a matrix keeps positive eigenvalues
  Ysample3 <- normalize_trace(Ysample)
  expect_equal(diag(mmean(Ysample3)), rep(1/3, 3), tol = 1E-2)
  expect_equal(eigen_desc(mmean(Ysample3))$values, rep(1/3, 3), tol = 1E-2)
  # bootstrapping large sample takes time:
  expect_gt(test_multiplicity(Ysample3, 3)$pval, 0.05)
  expect_gt(test_multiplicity_nonnegative(Ysample3, 3)$pval, 0.05)
})

rftiso <- function(n, meanlog = log(1/3) - 0.1^2/2, sdlog = 0.1){
  Ysample <- rsymm_lognorm(n, meanlog = diag(meanlog, 3), sigmalog = diag(sdlog^2, 6))
  normalise_trace(Ysample)
}

test_that("test on norm has bad size for small n, even with well confined sample eigenvectors and no normalization", {
  #  this test checks that poor convex hull behaviour is not due to non-negative distributions (I suspect it is due to eigenvectors with uniform distribution). Note that the transformation-based bootstrap does not have this problem.
  skip_if_fast_check() #takes a few seconds and not essential
  vals <- lapply(1:100, function(seed){
    set.seed(seed)
    Ysample <- rsymm_norm(15, mean = diag(1/3, 3), sigma = diag(c(0.001, 0.0001, 0.0001, 0.001, 0.0001, 0.001)))
    res <- suppressWarnings(test_multiplicity_nonnegative(Ysample, mult = 3, B = 10))
    c(res[c("pval", "t0", "B")], list(Ysample = Ysample))
  })
  # large statistic values seem associated to being outside the convex hull
  inhull <- !is.na(unlist(lapply(vals, "[[", 3)))
  expect_gt(mean(!inhull), 0.07) #about 10% usually
})

test_that("test has uniform distribution at large n=30", {
  set.seed(10)
  vals <- replicate(ifelse(fast_check_on(), 10, 100), {
    Ysample <- rftiso(30)
    res <- test_multiplicity_nonnegative(Ysample, mult = 3, B = 100)
    res[c("pval", "B")]
  })
  expect_equal(sum(is.na(unlist(vals[2, ]))), 0)
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings(ks.test(unlist(vals["pval",]), "punif"))
  expect_gt(res$p.value, 0.2)
})


test_that("test rejects some incorrect hypotheses p = 3", {
  skip_if_fast_check()
  set.seed(13321)
  Ysample <- rsymm(30, diag(c(2,2,1)-3), sigma = diag(3 * (3+1) /2))
  nonnegevals <- t(apply(Ysample, 1, function(v){
    es <- eigen_desc(invvech(v))
    vech(es$vectors %*% diag(exp(es$values)) %*% t(es$vectors))
  }))
  nonnegevals <- as_fsm(nonnegevals)
  eigen_desc(mmean(nonnegevals))$values
  set.seed(3654)
  res <- test_multiplicity_nonnegative(nonnegevals, mult = c(2,1), B = 1000)
  expect_gt(res$pval, 0.1)
  
  #test power is poor for the following situation
  #set.seed(3543)
  #expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(1, 2), B = 1000)$pval, 0.05)
  
  set.seed(3541)
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3), B = 1000)$pval, 0.05)
})



test_that("fixed trace from projection preserved by standardisation and ignored by stat", {
  set.seed(134)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(1000, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2))
  Ysample_n <- project_trace(Ysample)
  
  # check that only first element of Hevals  is changed by trace fix
  evals <- eigen_desc(mmean(Ysample))$values
  evals_n <- eigen_desc(mmean(Ysample_n))$values
  expect_equal((helmert(sum(mult)) %*% evals)[-1],
    (helmert(sum(mult)) %*% evals_n)[-1])
  
  std <- standardise_multiplicity(Ysample, mult)
  std_n <- standardise_multiplicity(Ysample_n, mult)
  expect_true(all.equal(as_fsm(apply(std, 1, function(m){projtrace_matrix(invvech(m))}, simplify = FALSE)), 
                        std_n, check.attributes = FALSE))
  
  set.seed(123)
  stat_orig <- stat_multiplicity(Ysample, mult = mult)
  set.seed(123)
  stat_n <- stat_multiplicity(Ysample_n, mult = mult)
  expect_equal(stat_orig, stat_n)
})
