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

test_that("tauhat and scalehat get close really to correct tau", {
  s = 2
  tau = 1/4
  p = 3
  covmat <- covOI(p, s, tau)
  set.seed(344)
  ms <- rsymm_norm(1E5, mean = diag(c(4,2,1)), sigma = covmat)
  Mhat <- invvech(colMeans(ms))
  tauest <- tauhat(ms, Mhat)
  
  expect_equal(attr(tauest, "numerator"),
    (1-p*(p+1)/2) * p * (tau/(1-p*tau)) * s^2,
    tolerance = 1E-2)
  expect_equal(attr(tauest, "denominator"),
               (p*(p+1)/2 - 1) * p * (1 + p * tau/(1-p*tau)) * s^2,
               tolerance = 1E-2)
  expect_equal(as.numeric(tauest), tau, tolerance = 1E-3, ignore_attr = TRUE)
  
  ssq <- scalesqhat(ms = ms, Mhat = Mhat)
  expect_equal(ssq, s^2, tolerance = 1E-2)
})
  