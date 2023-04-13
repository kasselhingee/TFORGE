test_that("vecd() works", {
  expect_equal(vecd(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1, 5, 9, sqrt(2) * c(2,3,6)))
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
  p <- 3
  C2 <- C1 <- diag(p*(p+1)/2)
  
  # set up distribution means
  M0 <- diag(c(1,2,4))
 
  Sigma <- rWishart(1, 6, diag(6)) 
  # simulate test statistic
  simulateTstat <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- stat_schwartzmann_eval(ms1, ms2)
    return(res$t)
  }
  simts <- replicate(100, simulateTstat(50, 100))

  #map from vech back to a matrix and then forward with vecd
  #vecd(invvech(1:6))
  #vecd(invvech(diag(C1)))
  
  # approximate distribution parameters
  anv <- S_anv(50, 100, mn1, mn2, 
               C1 = diag(vecd(invvech(diag(C1)))), 
               C2 = diag(vecd(invvech(diag(C1))))) #the extra work accounts for strange vecd of Schwartzmann
  qqplot(simts/anv$a,
    qchisq(ppoints(100), df = anv$v))
})

test_that("S_anv() gives exact distribution for 'OI' covariances", {
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
    return(res$p)
  }
  simts <- replicate(1000, simulateTstat(50, 100))
  plot(ecdf(simts))  #pvalue distribution should be approximately uniform, but it isn't
})
