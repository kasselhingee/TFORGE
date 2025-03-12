test_that("rsymm_lognorm normalised has correct eigenvalues", {
  skip_on_cran() #very slow to test
  # simulating from a narrow distribution to keep within bounds
  sdlog <- 0.1
  meanlogval <- log(1/3) - sdlog^2/2
  set.seed(3)
  Ysample <- rsymm_lognorm(1E5, meanlog = diag(meanlogval, 3), sigmalog = diag(sdlog^2, 6))
  
  #without projection/normalisation, the mean is slightly too high
  expect_true(all(eigen_desc(mmean(Ysample))$values > 1/3))

  # projection of matrix doesnt neccesarily keep positive eigenvalues?!
  Ysample2 <- project_trace(Ysample)
  Ysample2[, c(1, 4, 6)] <- Ysample2[, c(1, 4, 6)] + 1/3
  expect_equal(diag(mmean(Ysample2)), rep(1/3, 3), tol = 1E-3)
  expect_gt(suppressWarnings(test_fixedtrace(Ysample2, evals = rep(1/3, 3), B = 1000)$pval), 0.05)
  expect_gt(test_multiplicity(Ysample2, 3)$pval, 0.05)
  expect_gt(test_multiplicity_nonnegative(Ysample2, 3)$pval, 0.05)
  
  # normalisation of a matrix does though!
  Ysample3 <- normalize_trace(Ysample)
  expect_equal(diag(mmean(Ysample3)), rep(1/3, 3), tol = 1E-3)
  expect_equal(eigen_desc(mmean(Ysample3))$values, rep(1/3, 3), tol = 1E-3)
  expect_gt(test_multiplicity(Ysample3, 3)$pval, 0.05)
  expect_gt(test_multiplicity_nonnegative(Ysample3, 3)$pval, 0.05)
})

rftiso <- function(n, meanlog = log(1/3) - 0.1^2/2, sdlog = 0.1){
  Ysample <- rsymm_lognorm(n, meanlog = diag(meanlog, 3), sigmalog = diag(sdlog^2, 6))
  normalise_trace(Ysample)
}


test_that("test on lognorm has bad size for small n", {
  set.seed(10)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rftiso(15)
    res <- test_multiplicity_nonnegative(Ysample, mult = 3, B = 1000)
    res[c("pval", "B")]
  }, cl = 3)
  expect_gt(mean(is.na(unlist(vals[2, ]))), 0.07) # not in convex hull 0.09 of the time!
  pvals <- unlist(vals[1, ])
  expect_gt(abs(mean(pvals <= 0.05)-0.05), 0.01)
  qqplot(pvals, y = runif(1000))
  ks.test(pvals, "punif")
})

test_that("test on norm has bad size for small n", {
  set.seed(10)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rsymm_norm(15, mean = diag(1/3, 3), sigma = diag(0.1^2, 6))
    Ysample <- normalize_trace(Ysample)
    res <- test_multiplicity_nonnegative(Ysample, mult = 3, B = 1000)
    res[c("pval", "B")]
  }, cl = 3)
  expect_gt(mean(is.na(unlist(vals[2, ]))), 0.07) # not in convex hull 0.09 of the time!
  pvals <- unlist(vals[1, ])
  expect_gt(abs(mean(pvals <= 0.05)-0.05), 0.01)
  qqplot(pvals, y = runif(1000))
  ks.test(pvals, "punif")
})

test_that("test on norm has bad size for small n, even with well confined sample eigenvectors and no normalization", {
  vals <- pbapply::pblapply(1:1000, function(seed){
    set.seed(seed)
    Ysample <- rsymm_norm(15, mean = diag(1/3, 3), sigma = diag(c(0.001, 0.0001, 0.0001, 0.001, 0.0001, 0.001)))
    res <- test_multiplicity_nonnegative(Ysample, mult = 3, B = 1000)
    c(res[c("pval", "t0", "B")], list(Ysample = Ysample))
  }, cl = 3)
  # large statistic values seem associated to being outside the convex hull - strange
  inhull <- !is.na(unlist(lapply(vals, "[[", 3)))
  expect_gt(mean(!inhull), 0.07) #about 10%
})

test_that("bad size unique to weighted bootstrapping", {
  set.seed(10)
  vals <- pbapply::pbreplicate(100, {
    Ysample <- rsymm_norm(15, mean = diag(1/3, 3), sigma = diag(0.1^2, 6))
    Ysample <- normalize_trace(Ysample)
    res <- test_multiplicity(Ysample, mult = 3, B = 1000)
    res[c("pval", "B")]
  }, cl = 3)
  pvals <- unlist(vals[1, ])
  expect_lt(abs(mean(pvals <= 0.05)-0.05), 0.01)
  qqplot(pvals, y = runif(1000))
  expect_gt(suppressWarnings(ks.test(pvals, "punif"))$p.value, 0.05)
})

test_that("test has uniform distribution at large n", {
  set.seed(10)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rftiso(100)
    res <- test_multiplicity_nonnegative(Ysample, mult = 3, B = 1000)
    res[c("pval", "B")]
  }, cl = 3)
  expect_equal(sum(is.na(unlist(vals[2, ]))), 0)
  pvals <- unlist(vals[1, ])
  expect_lt(abs(mean(pvals <= 0.05)-0.05), 0.01)
  qqplot(pvals, y = runif(1000))
  ks.test(pvals[1:100], "punif")
})
