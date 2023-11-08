test_that("covOI gives correct matrices", {
  expect_equal(covOI(3, 2, 0, vectorisor = "vecd"),
               4 * diag(6))
  expect_equal(covOI(3, 1, 1/4, vectorisor = "vecd"),
               blockdiag(1 + diag(3), diag(3)))
  
  vechcovOI <- covOI(3, 1, 1/4, vectorisor = "vech")
  expect_equal(diag(vechcovOI), c(2,0.5,0.5,2,0.5,2))
  expect_equal(which(vechcovOI == 1, arr.ind = TRUE, useNames = FALSE),
  matrix(c(4,1,
           6,1,
           1,4,
           6,4,
           1,6,
           4,6), ncol = 2, byrow = TRUE), ignore_attr = TRUE)
})

test_that("OIinnerprod fast matches slow method", {
  s = 2
  tau = 1/4
  p = 3
  A <- invvech(rsymm_norm(1, mean = diag(3))[1, ])[[1]]
  B <- invvech(rsymm_norm(1, mean = diag(3))[1, ])[[1]]
  
  covmat <- covOI(3, s, tau, vectorisor = "vecd")
  slowinnprod <- drop(vecd(A) %*% solve(covmat) %*% vecd(B))
  fastinnprod <- OIinnerprod(A, B, s, tau)
  expect_equal(slowinnprod, fastinnprod) 
  
  fastinnprod2 <- OIinnerprod_sst(as.sst(list(A)), as.sst(list(B)), s, tau)
  expect_equal(fastinnprod, fastinnprod2)
  
  fastinnprod_twice <- OIinnerprod_sst(as.sst(list(A, A)), as.sst(list(B, B)), s, tau)
  expect_equal(fastinnprod_twice, c(fastinnprod2, fastinnprod2))
})

test_that("estimateOIparams get close really to correct tau and scale", {
  s = 2
  tau = 1/4
  p = 3
  covmat <- covOI(p, s, tau)
  set.seed(344)
  ms <- rsymm_norm(1E5, mean = diag(c(4,2,1)), sigma = covmat)
  Mhat <- invvech(colMeans(ms))
  OIparams <- estimateOIparams(ms, Mhat)
  expect_equal(OIparams$tau, tau, tolerance = 1E-3, ignore_attr = TRUE)
  expect_equal(OIparams$scalesq, s^2, tolerance = 1E-2)
})
  
test_that("testOIcov has uniform p values for a null situation", {
  s = 1
  tau = 1/8
  p = 3
  covmat <- covOI(p, s, tau, vectorisor = "vech")
  set.seed(344)
  vals <- pbapply::pbreplicate(1E2,
    {
    ms <- rsymm_norm(1E4, mean = diag(c(4,2,1)), sigma = covmat)
    res <- testOIcov(ms)
    c(pval = res$pval, stat = res$stat)
    })
  # hist(vals["stat", ])
  qqplot(vals["pval", ], y = runif(1000))
  expect_gt(ks.test(vals["pval", ], "punif")$p.value, 0.2)

  # looks like stat is missing a 6! 
  qqplot(vals["stat", ] + 6, y = rchisq(1000, df = 6 * (6+1)/2 -2)); abline(a=0, b= 1)
  offsetpvals <-  1 - pchisq(vals["stat", ]+6, df = 6 * (6+1)/2 -2)
  qqplot(offsetpvals, y = runif(1000)); abline(a=0, b= 1)
  expect_gt(ks.test(offsetpvals, "punif")$p.value, 0.2)
})

test_that("blk() returns correct averages", {
  evals <- 10:1
  mult <- c(2,2,2,2,2)
  expect_equal(blk(evals, mult), c(9.5,9.5, 7.5,7.5, 5.5,5.5, 3.5,3.5, 1.5,1.5))
  expect_equal(blk(evals, c(5,5)), c(rep(8, 5), rep(3, 5)))
  expect_equal(blk(evals, c(5,4,1)), c(rep(8, 5), rep(14/4, 4), 1))
})

test_that("test_multiplicity_OI() on null has uniform p values", {
  set.seed(13313)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rsymm_norm(1E3, diag(evals), sigma = covOI(length(evals), 1/2, 0, vectorisor = "vech"))
    test_multiplicity_OI(Ysample, mult = mult)$pval
  })
  
  qqplot(vals, y = runif(1000))
  res <- ks.test(vals, "punif")
  expect_gt(res$p.value, 0.15)
})
  