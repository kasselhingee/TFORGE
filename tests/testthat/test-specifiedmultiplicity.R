test_that("stat_specifiedmultiplicity() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  av <- mmean(Ysample)
  es <- eigen(av)
  Ystdsample <- standardise_specifiedmultiplicity(Ysample, mult = c(2, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen(newav)
  expect_equal(newes$values[1], newes$values[2])
  expect_equal(newes$values[1], 2.5, tolerance = 0.3)
  expect_equal(newes$values[3], es$values[3])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(2,1)), 0)
  expect_error(expect_equal(stat_specifiedmultiplicity(Ysample, mult = c(2,1)), 0))
})

test_that("stat_specifiedmultiplicity() is zero for standarised sample, dim 7", {
  set.seed(13131)
  Ysample <- rsymm(50, 7)
  Ysample <- lapply(Ysample, `+`, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  av <- mmean(Ysample)
  es <- eigen(av)
  Ystdsample <- standardise_specifiedmultiplicity(Ysample, mult = c(3, 2, 1, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen(newav)
  expect_equal(newes$values[1:3], rep(3, 3), tolerance = 0.1)
  expect_equal(newes$values[4:5], rep(2, 2), tolerance = 0.1)
  expect_equal(newes$values[6], es$values[6])
  expect_equal(newes$values[7], es$values[7])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(3, 2, 1, 1)), 0)
  expect_error(expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(3, 2, 2)), 0))
  expect_error(expect_equal(stat_specifiedmultiplicity(Ysample, mult = c(3, 2, 1, 1)), 0))
})

test_that("stat_specifiedevals() doesn't reject for simulation of single sample from null, and rejects otherwise", {
  set.seed(1331)
  Ysample <- rsymm(50, 7)
  Ysample <- lapply(Ysample, `+`, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  res <- test_specifiedmultiplicity(Ysample, mult = c(3,2,1,1), 100)
  expect_gt(res$pval, 0.2)
  
  res <- test_specifiedmultiplicity(Ysample, mult = c(3,2,2), 100)
  expect_lt(res$pval, 0.2)
})

test_that("xiget() behaves", {
  p = 5
  n = 10
  mult <- c(3, 2)
  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  
  expect_equal(xiget(c(3.0, 3.0, 3, 2, 2), mult, idxs), rep(0, 3))
  expect_equal(xiget(c(3.0, 3.0, 3.1, 2, 2), mult, idxs)==0, c(TRUE, FALSE, TRUE))
  expect_equal(xiget(c(3.0, 3.0, 3.1, 2, 2.1), mult, idxs)==0, c(TRUE, FALSE, FALSE))
  expect_equal(xiget(c(3.0, 2.9, 3.1, 2, 2.1), mult, idxs)==0, c(FALSE, FALSE, FALSE))
})

test_that("xicovar() gives the same covariance as sample covariance", {
  p = 5
  n = 10
  mult <- c(3, 2)
  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  
  simxi <- function(n, mn, sigma, mult, idxs, evecs = NULL){
    Ysample <- rsymm(n, mn, sigma)
    Ybar <- mmean(Ysample)
    if (is.null(evecs)){evals <- eigen(Ybar)$values} 
    else {evals <- diag(t(evecs) %*% Ybar %*% evecs)} #use the true eigenvectors
    xi <- xiget(evals, mult, idxs)
    return(xi)
  }
  
  # set up distribution covariance
  set.seed(316)
  nvar <- (p + 1) * p/2
  C0_U <- mclust::randomOrthogonalMatrix(nvar, nvar)
  C0 <- C0_U %*% diag(1:nrow(C0_U)) %*% t(C0_U) #use this to simulate
  set.seed(348)
  mn_U <- mclust::randomOrthogonalMatrix(p, p)
  mn <- mn_U %*% diag(c(3, 3, 3, 2, 2)) %*% t(mn_U)
  
  # theoretical covariance
  thecov <- xicovar(mult, idxs, eigen(mn)$vectors, C0/n)
 
  set.seed(35468) 
  emcov <- replicate(100,
   simxi(n, mn = mn, sigma = C0, mult, idxs, eigen(mn)$vectors)) |>
    t() |>
    cov()
  expect_equal(emcov, thecov)
  abs(emcov - thecov) / abs(thecov)
})
