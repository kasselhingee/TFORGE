test_that("vecd() works", {
  expect_equal(vecd(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1, 5, 9, sqrt(2) * c(2,3,6)))
  
  m <- rsymm(1, diag(3))[[1]]
  expect_equal(invvecd(vecd(m)), m)
})

test_that("S_mcovar gets close to true", {
  skip("slow")
  Ysample <- rsymm(100000, matrix(0, 3, 3), sigma = diag( 3 * (3+1)/2))
  expect_silent(c0 <- S_mcovar(merr(Ysample)))
  expect_equal(c0, diag(c(rep(1, 3), rep(2, 3))), tolerance = 0.06)
})

test_that("makeblockdiagonal works", {
  A <- matrix(1, 3, 5)
  B <- matrix(2, 4, 6)
  expect_equal(blockdiag(A, B)[4,1], 0)
  expect_equal(blockdiag(A, B)[4,6], 2)
})

test_that("Results are consistent with the simulation check at the start of Schwartzmann Section 3.1.", {
  # set up distribution means
  M0 <- diag(c(1,2,4))
 
  simulatestat <- function(n1, n2, M0){
    # simulate Sigma
    Sigma <- drop(rWishart(1, 6, diag(6)))
    #simulate samples
    ms1 <- rsymm_Schwartzmann(n1, M0, Sigma)
    ms2 <- rsymm_Schwartzmann(n2, M0, Sigma)
    res <- stat_schwartzmann_eval(ms1, ms2)
    return(unlist(res))
  }
  
  # fast
  set.seed(1354)
  sims <- replicate(100, simulatestat(50, 50, M0))
  ecdffun <- ecdf(sims["pval", ])
  #plot(ecdffun); abline(0, 1, lty = "dotted") #pvalue distribution should be approximately uniform for the null
  expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.1)
  
  # much slower
  # set.seed(3546)
  # sims <- replicate(10000, simulatestat(50, 50, M0))
  # ecdffun <- ecdf(sims["pval", ])
  # plot(ecdffun); abline(0, 1, lty = "dotted")
  # expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.05)
})

test_that("Results are consistent with the simulation check at the start of Schwartzmann Section 3.1 with fixed Sigma", {
  # set up distribution means
  M0 <- diag(c(1,2,4))

  # simulate Sigma  
  set.seed(3145)
  Sigma <- drop(rWishart(1, 6, diag(6)))
  
  simulatestat <- function(n1, n2, M0){
    #simulate samples
    ms1 <- rsymm_Schwartzmann(n1, M0, Sigma)
    ms2 <- rsymm_Schwartzmann(n2, M0, Sigma)
    res <- stat_schwartzmann_eval(ms1, ms2)
    return(unlist(res))
  }
  
  # fast
  set.seed(1354)
  sims <- replicate(100, simulatestat(50, 50, M0))
  ecdffun <- ecdf(sims["pval", ])
  #plot(ecdffun); abline(0, 1, lty = "dotted") #pvalue distribution should be approximately uniform for the null
  expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.1)
  
  # much slower
  # set.seed(3546)
  # sims <- replicate(10000, simulatestat(50, 50, M0))
  # ecdffun <- ecdf(sims["pval", ])
  # plot(ecdffun); abline(0, 1, lty = "dotted")
  # expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.05)
})

test_that("S_anv() gives good distribution of Tstat for diagonal covariances", {
  skip("Schwartzmann's distribution approximation using two moments and projection onto a manifold doesn't appear to converge with increasing sample sizes.")
  p <- 3
  C2 <- C1 <- diag(p*(p+1)/2)
  
  # set up distribution means
  set.seed(348)
  mn_U1 <- mclust::randomOrthogonalMatrix(p, p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(543)
  mn_U2 <- mclust::randomOrthogonalMatrix(p, p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # simulate test statistic
  simulateTstat <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- stat_schwartzmann_eval(ms1, ms2)
    return(unlist(res))
  }
  set.seed(231654)
  sims <- replicate(100, simulateTstat(500, 1000))
  # ecdffun <- ecdf(sims["pval", ])
  # plot(ecdffun); abline(0, 1, lty = "dotted") #pvalue distribution should be approximately uniform for the null
  expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.1)
  
  # distribution parameters
  anv <- S_anv(500, 1000, mn1, mn2, 
               C1 = diag(vecd(invvech(diag(C1)))), 
               C2 = diag(vecd(invvech(diag(C1))))) #the extra work accounts for strange vecd of Schwartzmann
  
  qqplot(sims["t", ]/anv$a,
         qchisq(ppoints(100), df = anv$v))
  abline(0, 1, lty = "dotted")
  
  plot(ecdf(pchisq(sims["t", ]/anv$a, df = anv$v)))
  abline(0, 1, lty = "dotted")
})

test_that("Schwartzmann's Statistic (Tstatstar) is well approximated by S_anv() with a chisq", {
  skip("Schwartzmann's distribution approximation using two moments doesn't appear to converge with increasing sample sizes.")
  p <- 3
  C2 <- C1 <- diag(p*(p+1)/2)
  n1 <- 1000 #Schwartzmann's distribution approximation using two moments doesn't appear to converge with increasing sample sizes.
  n2 <- 1000
  # set up distribution means
  set.seed(1354)
  mn_U1 <- mclust::randomOrthogonalMatrix(p, p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(135)
  mn_U2 <- mclust::randomOrthogonalMatrix(p, p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # distribution parameters
  anv <- S_anv(n1, n2, mn1, mn2, 
               C1 = diag(vecd(invvech(diag(C1)))), 
               C2 = diag(vecd(invvech(diag(C1)))))
  
  # more exact statistic simulation (still only an approximation of the distribution!)
  simulateTstatstar <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- statstar_schwartzmann_eval(ms1, ms2, mn1, mn2)
    return(unlist(res))
    }
  set.seed(31456)
  sims <- replicate(1000, simulateTstatstar(n1, n2))
  
  qqplot(sims/anv$a,
         qchisq(ppoints(1000), df = anv$v))
  abline(0, 1, lty = "dotted")
  plot(ecdf(pchisq(sims/anv$a, df = anv$v))); abline(0, 1, lty = "dotted")
})


test_that("stat_schwartzmann_eval() doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(3,2,1)))
  )
  
  res <- stat_schwartzmann_eval(Ysamples[[1]], Ysamples[[2]])
  expect_gt(res$pval, 0.2)
})

test_that("stat_schwartzmann_eval() reject for simulation of multi sample not from null", {
  set.seed(13)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(4,3,2)))
  )
  
  res <- stat_schwartzmann_eval(Ysamples[[1]], Ysamples[[2]])
  expect_lt(res$pval, 0.05)
})
