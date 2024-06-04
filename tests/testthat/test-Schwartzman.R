test_that("vecd() works", {
  expect_equal(vecd(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1, 5, 9, sqrt(2) * c(2,3,6)))
  
  m <- invvech(rsymm(1, diag(3))[1, ])
  expect_equal(invvecd(vecd(m)), m)
})
test_that("vech2vecd works", {
  set.seed(354)
  m <- invvech(rsymm(1, diag(3))[1, ])
  expect_equal(invvecd(vech2vecd(vech(m))), m)
})

test_that("S_mcovar gets close to true", {
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

test_that("Results are consistent with the simulation check at the start of Schwartzman Section 3.1.", {
  # set up distribution means
  M0 <- diag(c(1,2,4))
 
  simulatestat <- function(n1, n2, M0){
    # simulate Sigma
    Sigma <- drop(rWishart(1, 6, diag(6)))
    #simulate samples
    ms1 <- rsymm_Schwartzman(n1, M0, Sigma)
    ms2 <- rsymm_Schwartzman(n2, M0, Sigma)
    res <- stat_schwartzman_eval(ms1, ms2)
    return(unlist(res))
  }
  
  # fast
  set.seed(1354)
  sims <- replicate(100, simulatestat(50, 50, M0))
  ecdffun <- ecdf(sims["pval", ])
  #plot(ecdffun); abline(0, 1, lty = "dotted") #pvalue distribution should be approximately uniform for the null
  expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.1)
  expect_gt(ks.test(sims["pval", ], "punif")$p.value, 0.2)
  
  # much slower
  # set.seed(3546)
  # sims <- replicate(10000, simulatestat(50, 50, M0))
  # ecdffun <- ecdf(sims["pval", ])
  # plot(ecdffun); abline(0, 1, lty = "dotted")
  # expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.05)
})

test_that("Results are consistent with the simulation check at the start of Schwartzman Section 3.1 with fixed Sigma", {
  # set up distribution means
  M0 <- diag(c(1,2,4))

  # simulate Sigma  
  set.seed(3145)
  Sigma <- drop(rWishart(1, 6, diag(6)))
  
  simulatestat <- function(n1, n2, M0){
    #simulate samples
    ms1 <- rsymm_Schwartzman(n1, M0, Sigma)
    ms2 <- rsymm_Schwartzman(n2, M0, Sigma)
    res <- stat_schwartzman_eval(ms1, ms2)
    return(unlist(res))
  }
  
  # fast
  set.seed(1354)
  sims <- replicate(100, simulatestat(50, 50, M0))
  ecdffun <- ecdf(sims["pval", ])
  #plot(ecdffun); abline(0, 1, lty = "dotted") #pvalue distribution should be approximately uniform for the null
  expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.1)
  expect_gt(ks.test(sims["pval", ], "punif")$p.value, 0.2)
  
  # much slower
  # set.seed(3546)
  # sims <- replicate(10000, simulatestat(50, 50, M0))
  # ecdffun <- ecdf(sims["pval", ])
  # plot(ecdffun); abline(0, 1, lty = "dotted")
  # expect_lt(max(abs(ecdffun(ppoints(100)) - ppoints(100))), 0.05)
})

test_that("Schwartzman statistic for iid elements is chisq under H0 and pvals uniform", {
  set.seed(15)
  vals <- replicate(1000, {
    Ysamples <- list(
      rsymm(100, diag(c(3,2,1))), #50 was too small, 100 is needed to get good estimates of a and v for the chisquared distribution
      rsymm(100, diag(c(3,2,1)))
    )
    res <- stat_schwartzman_eval(Ysamples[[1]], Ysamples[[2]])
    res})

  # qqplot(unlist(vals["t", ]), y = rchisq(1000, df = 3)); abline(0, 1, lty = "dotted")
  res <- ks.test(unlist(vals["t", ]), "pchisq", df = 3)
  expect_gt(res$p.value, 0.2)
  
  # true value of a is 1
  # true value of v is 3
  anv <- S_anv(100,100,diag(c(3,2,1)),diag(c(3,2,1)),
               C1 = vech2vecd_mat(6) %*% diag(6) %*% t(vech2vecd_mat(6)),
               C2 = vech2vecd_mat(6) %*% diag(6) %*% t(vech2vecd_mat(6)))
  expect_equal(anv$a, 1)
  expect_equal(anv$v, 3)
  
  # hist(unlist(vals["a", ])) #biased - tends to be higher than 1
  expect_equal(mean(unlist(vals["a", ])), 1, tolerance = 0.1)
  # hist(unlist(vals["v", ])) #biased - always less than 3
  expect_equal(mean(unlist(vals["v", ])), 3, tolerance = 0.1)
  expect_equal(max(unlist(vals["v", ])), 3, tolerance = 1E-2)
  
  qqplot(unlist(vals["pval", ]), y = runif(1000)); abline(0, 1, lty = "dotted") #it looks very uniform!
  res <- ks.test(unlist(vals["pval", ]), "punif")
  expect_gt(res$p.value, 0.2)
})


test_that("pvalue close to uniform for non diagonal mean, unequal sample size", {
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
    res <- stat_schwartzman_eval(ms1, ms2)
    return(unlist(res))
  }
  set.seed(231654)
  sims <- replicate(1000, simulateTstat(100, 500))
  # qqplot(sims["pval", ], runif(1000)); abline(0, 1, lty = "dotted")
  expect_gt(ks.test(sims["pval", ], "punif")$p.value, 0.2)
})

test_that("pvalue close to uniform when some correlation", {
  p <- 3
  C2 <- C1 <- matrix(0.3, 6, 6)
  diag(C2) <- diag(C1) <- 1
  
  # set up distribution means
  set.seed(3481)
  mn_U1 <- mclust::randomOrthogonalMatrix(p, p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(5431)
  mn_U2 <- mclust::randomOrthogonalMatrix(p, p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # simulate test statistic
  simulateTstat <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- stat_schwartzman_eval(ms1, ms2)
    return(unlist(res))
  }
  set.seed(22)
  sims <- replicate(1000, simulateTstat(100, 500))
  # qqplot(sims["pval", ], runif(1000)); abline(0, 1, lty = "dotted")
  expect_gt(ks.test(sims["pval", ], "punif")$p.value, 0.2)
})


test_that("convergence of S_anv() with increasing sample size", {
  p <- 3
  C2 <- C1 <- matrix(0.6, 6, 6)
  diag(C2) <- diag(C1) <- 1
  # set up distribution means
  set.seed(34811)
  mn_U1 <- mclust::randomOrthogonalMatrix(p, p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(54311)
  mn_U2 <- mclust::randomOrthogonalMatrix(p, p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # population parameters factored out
  anv <- S_anv(1E6, 1E6, mn1, mn2, 
               C1 = vech2vecd_mat(6) %*% C1 %*% t(vech2vecd_mat(6)),
               C2 = vech2vecd_mat(6) %*% C2 %*% t(vech2vecd_mat(6)))
  
  # Sample covariance is Wishart distributed
  set.seed(354)
  ests <- lapply(c(10,100,1000), function(n1){
    n2 <- n1
    unlist(S_anv(n1, n2, mn1, mn2, 
          C1 = vech2vecd_mat(6) %*% drop(rWishart(1, n1, C1)/(n1-1)) %*% t(vech2vecd_mat(6)),
          C2 = vech2vecd_mat(6) %*% drop(rWishart(1, n2, C2)/(n2-1)) %*% t(vech2vecd_mat(6))))
  })
  estanv <- simplify2array(c(list(unlist(anv)), ests))
  err <- abs(estanv[,-1] - estanv[, 1])
  expect_equal(order(err["a", ], decreasing = TRUE), 1:3)
  expect_equal(order(err["v", ], decreasing = TRUE), 1:3)
  expect_lt(err["a", 3], 1E-1)
  expect_lt(err["v", 3], 1E-1)
})



test_that("(Welch-Satterthwaite approximation) Schwartzman's Tstatstar has moments implied by S_anv(), cov = I", {
  p <- 3
  C2 <- C1 <- solve(vech2vecd_mat(6)) %*% diag(6) %*% t(solve(vech2vecd_mat(6)))
  n1 <- 10 
  n2 <- 10
  # set up distribution means
  set.seed(1354)
  mn_U1 <- runifortho(p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(135)
  mn_U2 <- runifortho(p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # distribution parameters
  anv <- S_anv(n1, n2, mn1, mn2, 
               C1 = vech2vecd_mat(6) %*% C1 %*% t(vech2vecd_mat(6)),
               C2 = vech2vecd_mat(6) %*% C2 %*% t(vech2vecd_mat(6)))
  # when covariances C1=C2 are I and n1=n2=n then we have mean = tr(Lambda) is the below
  target_trLambda <- sum(diag(Omega_eval(n1, n2, mn_U1, mn_U2)))/n1
  expect_equal(anv$a * anv$v, target_trLambda)
  
  # simulate statstar, result should have mean and var given by anv
  simulateTstatstar <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- statstar_schwartzman_eval(ms1, ms2, mn1, mn2)
    return(unlist(res))
    }
  set.seed(31456)
  sims <- replicate(1E4, simulateTstatstar(n1, n2))
  # from Casella and Berger (Statistical Inference) on Satterthwaite's approximation
  # mean of statstar is sum(lamba_i * 1) = tr(Lambda). statstar/tr(Lambda) has mean of sum(lambda_i/tr(Lambda) * 1) = 1.
  # Satterthwaite estimated nu as:
  # hatnu = (sum( lambda_i/tr(Lambda) * Y_i)^2)/sum{ (lambda_i/tr(Lambda))^2 Y_i^2}
  # where Y_i = z_i^2.
  # Since we have E(Y_i) = 1 we can back pedal a bit and hatnu is actually
  # hatnu = 1 * tr(Lambda)^2 / sum{lambda_i^2} = tr(Lambda)^2/tr(Lambda^2)
  # that means statstar/tr(Lambda) d~d chisq(hatnu)/hatnu
  # so statstar d~d tr(Lambda)/hatnu * chisq(hatnu)/hatnu
  # so a = tr(Lambda) * tr(Lambda^2)/tr(Lambda)^2 = tr(Lambda^2)/tr(Lambda). This is exactly what Schwartzman writes.
  
  # the mean of sims should be tr(Lambda) * hatnu = anv$a/anv$nu
  expect_equal(mean(sims), anv$a * anv$v, tolerance = 2 * sd(sims)/sqrt(length(sims)))
  #similarly for the variance
  expect_equal(var(sims), 2*anv$a^2 * anv$v, tolerance = 0.1)
  
  # chisq approximation looks pretty good but definitely not exact
  # qqplot(sims/anv$a,
  #        qchisq(ppoints(1000), df = anv$v))
  # abline(0, 1, lty = "dotted")
})

test_that("(Welch-Satterthwaite approximation) Schwartzman's Tstatstar has moments implied by S_anv() general cov", {
  p <- 3
  set.seed(10)
  C1 <- drop(rWishart(1, 6, diag(6)))
  C2 <- drop(rWishart(1, 6, diag(6)))
  n1 <- 10 
  n2 <- 10
  # set up distribution means
  set.seed(1354)
  mn_U1 <- runifortho(p)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(135)
  mn_U2 <- runifortho(p)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # distribution parameters
  anv <- S_anv(n1, n2, mn1, mn2, 
               C1 = vech2vecd_mat(6) %*% C1 %*% t(vech2vecd_mat(6)),
               C2 = vech2vecd_mat(6) %*% C2 %*% t(vech2vecd_mat(6)))
  
  # simulate statstar, result should have mean and var given by anv
  simulateTstatstar <- function(n1, n2){
    ms1 <- rsymm(n1, mn1, C1)
    ms2 <- rsymm(n2, mn2, C2)
    res <- statstar_schwartzman_eval(ms1, ms2, mn1, mn2)
    return(unlist(res))
  }
  set.seed(31456)
  sims <- replicate(1E4, simulateTstatstar(n1, n2))
  # from Casella and Berger (Statistical Inference) on Satterthwaite's approximation
  # mean of statstar is sum(lamba_i * 1) = tr(Lambda). statstar/tr(Lambda) has mean of sum(lambda_i/tr(Lambda) * 1) = 1.
  # Satterthwaite estimated nu as:
  # hatnu = (sum( lambda_i/tr(Lambda) * Y_i)^2)/sum{ (lambda_i/tr(Lambda))^2 Y_i^2}
  # where Y_i = z_i^2.
  # Since we have E(Y_i) = 1 we can back pedal a bit and hatnu is actually
  # hatnu = 1 * tr(Lambda)^2 / sum{lambda_i^2} = tr(Lambda)^2/tr(Lambda^2)
  # that means statstar/tr(Lambda) d~d chisq(hatnu)/hatnu
  # so statstar d~d tr(Lambda)/hatnu * chisq(hatnu)/hatnu
  # so a = tr(Lambda) * tr(Lambda^2)/tr(Lambda)^2 = tr(Lambda^2)/tr(Lambda). This is exactly what Schwartzman writes.
  
  # the mean of sims should be tr(Lambda) * hatnu = anv$a/anv$nu
  expect_equal(mean(sims), anv$a * anv$v, tolerance = 2 * sd(sims)/sqrt(length(sims)))
  #similarly for the variance
  expect_equal(var(sims), 2*anv$a^2 * anv$v, tolerance = 0.1)
  
  # chisq approximation looks pretty good but definitely not exact
  # qqplot(sims/anv$a, qchisq(ppoints(1000), df = anv$v)); abline(0, 1, lty = "dotted")
})


test_that("Schwartzman's Omega(M) is correct for U1=U2=I", {
  out <- Omega_eval(1,1,diag(3), diag(3))
  #c(1,0,0,0,0,0, -1,0,0,0,0,0) %*% t(c(1,0,0,0,0,0, -1,0,0,0,0,0)) = 
  m <- diag(vecd(diag(3)))
  target <- 1/2 * rbind(cbind(m, -m),
                  cbind(-m, m))
  expect_equal(out, target)
})

test_that("Omega_eval returns correct trace", {
  set.seed(1354)
  mn_U1 <- runifortho(3)
  mn1 <- mn_U1 %*% diag(c(3,2,1)) %*% t(mn_U1)
  set.seed(135)
  mn_U2 <- runifortho(3)
  mn2 <- mn_U2 %*% diag(c(3,2,1)) %*% t(mn_U2)
  
  # trace of Omega is just sum of wi^2 for each i
  # wi is just halfish elements of U1[,i] %*% t(U1[,i]) with an extra factor of sqrt(2) for the off diagonals
  # so sum of just halfish elements of (U1[,i] %*% t(U1[,i]))^2 with an extra factor of 2 for the off diagonals
  target <- sum(vapply(1:3, function(i){
    m1 <- mn_U1[,i]^2 %*% t(mn_U1[,i]^2)
    m2 <- mn_U2[,i]^2 %*% t(mn_U2[,i]^2)
    m1[upper.tri(m1)] <- 2 * m1[upper.tri(m1)]
    m1[lower.tri(m1)] <- 0
    m2[upper.tri(m2)] <- 2 * m2[upper.tri(m2)]
    m2[lower.tri(m2)] <- 0
    sum(m1) + sum(m2)
  }, FUN.VALUE = 1.0))/2 #the 2 here is because 1*1/(1+1) = 2
  expect_equal(sum(diag(Omega_eval(1, 1, mn_U1, mn_U2))), target)
})

test_that("stat_schwartzman_eval() is invariant to rotations", {
  set.seed(15)
  Ysample1_0 <- rsymm(100, diag(c(0,0,0)), sigma = diag(6)/10)
  Ysample2 <- rsymm(100, diag(c(3,2,1)))
  Ysample1 <- t(t(Ysample1_0) + vech(diag(c(3,2,1))))
  set.seed(134)
  rotmat <- runifortho(3)
  Ysample1_r <- t(t(Ysample1_0) + vech(rotmat %*% diag(c(3,2,1)) %*% t(rotmat)))

  res1 <- stat_schwartzman_eval(Ysample1, Ysample2)
  res2 <- stat_schwartzman_eval(Ysample1_r, Ysample2)
  expect_equal(res1$pval, res2$pval, tolerance = 0.1)
})

test_that("stat_schwartzman_eval() doesn't reject for simulation of multi sample from null", {
  set.seed(15)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(3,2,1)))
  )
  
  res <- stat_schwartzman_eval(Ysamples[[1]], Ysamples[[2]])
  expect_gt(res$pval, 0.2)
  res2 <- stat_schwartzman_eval(as.mstorsst(Ysamples))
  expect_equal(res, res2)
})

test_that("stat_schwartzman_eval() reject for simulation of multi sample not from null", {
  set.seed(13)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(4,3,2)))
  )
  
  res <- stat_schwartzman_eval(Ysamples[[1]], Ysamples[[2]])
  expect_lt(res$pval, 0.05)
})
