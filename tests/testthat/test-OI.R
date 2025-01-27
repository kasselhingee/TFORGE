test_that("OIcov gives correct matrices", {
  expect_equal(OIcov(3, 2, 0, vectorisor = "vecd"),
               4 * diag(6))
  expect_equal(OIcov(3, 1, 1/4, vectorisor = "vecd"),
               blockdiag(1 + diag(3), diag(3)))
  
  vechOIcov <- OIcov(3, 1, 1/4, vectorisor = "vech")
  expect_equal(diag(vechOIcov), c(2,0.5,0.5,2,0.5,2))
  expect_equal(which(vechOIcov == 1, arr.ind = TRUE, useNames = FALSE),
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
  A <- invvech(rsymm_norm(1, mean = diag(3))[1, ])
  B <- invvech(rsymm_norm(1, mean = diag(3))[1, ])
  
  covmat <- OIcov(3, s, tau, vectorisor = "vecd")
  slowinnprod <- drop(vecd(A) %*% solve(covmat) %*% vecd(B))
  fastinnprod <- OIinnerprod(A, B, s, tau)
  expect_equal(slowinnprod, fastinnprod) 
  
  fastinnprod2 <- OIinnerprod_fsm(as_fsm(list(A)), as_fsm(list(B)), s, tau)
  expect_equal(fastinnprod, fastinnprod2)
  
  fastinnprod_twice <- OIinnerprod_fsm(as_fsm(list(A, A)), as_fsm(list(B, B)), s, tau)
  expect_equal(fastinnprod_twice, c(fastinnprod2, fastinnprod2))
})

test_that("estimate_OIcov get close really to correct tau and scale", {
  s = 2
  tau = 1/4
  p = 3
  covmat <- OIcov(p, s, tau)
  set.seed(344)
  ms <- rsymm_norm(1E5, mean = diag(c(4,2,1)), sigma = covmat)
  Mhat <- invvech(colMeans(ms))
  OIparams <- estimate_OIcov(ms, Mhat)
  expect_equal(OIparams$tau, tau, tolerance = 1E-3, ignore_attr = TRUE)
  expect_equal(OIparams$scalesq, s^2, tolerance = 1E-2)
})
  
test_that("blk() returns correct averages", {
  evals <- 10:1
  mult <- c(2,2,2,2,2)
  expect_equal(blk(evals, mult), c(9.5,9.5, 7.5,7.5, 5.5,5.5, 3.5,3.5, 1.5,1.5))
  expect_equal(blk(evals, c(5,5)), c(rep(8, 5), rep(3, 5)))
  expect_equal(blk(evals, c(5,4,1)), c(rep(8, 5), rep(14/4, 4), 1))
})

test_that("test_multiplicity_OI() on null has uniform p values", {
  set.seed(4)
  abasis <- runifortho(7)
  set.seed(13312)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- replicate(1000, {
    Ysample <- rsymm_norm(1E2, diag(evals), sigma = OIcov(length(evals), 1/2, 0, vectorisor = "vech"))
    c(r = test_multiplicity_OI(Ysample, mult = mult, refbasis = "random")$pval,
      c = test_multiplicity_OI(Ysample, mult = mult, refbasis = diag(1, 7))$pval,
      a = test_multiplicity_OI(Ysample, mult = mult, refbasis = abasis)$pval)
  })
  
  # qqplot(vals["r", ], y = runif(1000))
  # qqplot(vals["c", ], y = runif(1000))
  # qqplot(vals["a", ], y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals["r", ], "punif"))$p.value, 0.15)
  expect_gt(suppressWarnings(ks.test(vals["c", ], "punif"))$p.value, 0.15) 
  expect_gt(suppressWarnings(ks.test(vals["a", ], "punif"))$p.value, 0.15)
})

test_that("test_multiplicity_OI() bootstrap on null has uniform p values", {
  set.seed(4)
  abasis <- runifortho(7)
  set.seed(13312)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- pbapply::pbreplicate(100, {
    Ysample <- rsymm_norm(30, diag(evals), sigma = OIcov(length(evals), 1/2, 0, vectorisor = "vech"))
    c(r = test_multiplicity_OI(Ysample, mult = mult, B = 100, refbasis = "random")$pval,
      c = test_multiplicity_OI(Ysample, mult = mult, B = 100, refbasis = diag(1, 7))$pval,
      a = test_multiplicity_OI(Ysample, mult = mult, B = 100, refbasis = abasis)$pval)
  })
  
  # qqplot(vals["r", ], y = runif(1000))
  # qqplot(vals["c", ], y = runif(1000))
  # qqplot(vals["a", ], y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals["r", ], "punif"))$p.value, 0.15)
  expect_gt(suppressWarnings(ks.test(vals["c", ], "punif"))$p.value, 0.15) 
  expect_gt(suppressWarnings(ks.test(vals["a", ], "punif"))$p.value, 0.15)
})

