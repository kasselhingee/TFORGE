test_that("stat is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, mean = diag(c(3,2,1)))
  av <- mmean(Ysample)
  es <- eigen_desc(av)
  Ystdsample <- standardise_multiplicity(Ysample, mult = c(2, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen_desc(newav)
  expect_equal(newes$values[1], newes$values[2])
  expect_equal(newes$values[1], 2.5, tolerance = 0.3)
  expect_equal(newes$values[3], es$values[3])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_multiplicity(Ystdsample, mult = c(2,1)), 0)
  expect_error(expect_equal(stat_multiplicity(Ysample, mult = c(2,1)), 0))
})

test_that("stat is zero for standarised sample, dim 7", {
  set.seed(13131)
  Ysample <- rsymm(50, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  av <- mmean(Ysample)
  es <- eigen_desc(av)
  Ystdsample <- standardise_multiplicity(Ysample, mult = c(3, 2, 1, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen_desc(newav)
  expect_equal(newes$values[1:3], rep(3, 3), tolerance = 0.1)
  expect_equal(newes$values[4:5], rep(2, 2), tolerance = 0.1)
  expect_equal(newes$values[6], es$values[6])
  expect_equal(newes$values[7], es$values[7])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 1, 1)), 0)
  expect_error(expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 2)), 0))
  expect_error(expect_equal(stat_multiplicity(Ysample, mult = c(3, 2, 1, 1)), 0))
})

test_that("debugging stat with true evecs is same regardless of sortevecs", {
  set.seed(1332)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(30, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
  unsort <- suppressWarnings(stat_multiplicity(Ysample, mult = mult, evecs = diag(sum(mult)), sortevecs = FALSE))
  sort <- suppressWarnings(stat_multiplicity(Ysample, mult = mult, evecs = diag(sum(mult)), sortevecs = TRUE))
  
  expect_equal(sort, unsort)  
})

test_that("debugging stat with true evecs has correct null distribution", {
  set.seed(1332)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rsymm_norm(200, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    stat_multiplicity(Ysample, mult = mult, evecs = diag(sum(mult)))
  }, cl = 2)
  
  # qqplot(vals, y = rchisq(1000, df = sum(mult-1))) 
  res <- ks.test(vals, "pchisq", df = sum(mult-1))
  expect_gt(res$p.value, 0.2)
})

test_that("stat has correct null distribution", {
  set.seed(1331)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rsymm_norm(1000, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    stat_multiplicity(Ysample, mult = mult)
  })
  
  # qqplot(vals, y = rchisq(1000, df = sum(mult-1))) #looks more like df=5
  res <- ks.test(vals, "pchisq", df = sum(mult-1))
  expect_gt(res$p.value, 0.2)
})


test_that("test has uniform distribution", {
  skip_on_cran() #test very slow
  set.seed(1331)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- pbapply::pbreplicate(1000, {
    Ysample <- rsymm_norm(1000, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    test_multiplicity(Ysample, mult = mult, B = 100)$pval
  }, cl = 2)
  
  # qqplot(vals, y = runif(1000))
  res <- suppressWarnings(ks.test(vals, "punif"))
  expect_gt(res$p.value, 0.2)
})

test_that("test_multiplicity() doesn't reject for simulation of single sample from null, and rejects otherwise", {
  set.seed(1331)
  Ysample <- rsymm(100, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  res <- test_multiplicity(Ysample, mult = c(3,2,1,1), 100)
  expect_gt(res$pval, 0.2)
  
  expect_lt(test_multiplicity(Ysample, mult = c(3,3,1), 100)$pval, 0.05)
  expect_lt(test_multiplicity(Ysample, mult = c(2,3,1,1), 100)$pval, 0.05)
  expect_lt(test_multiplicity(Ysample, mult = c(2,3,2), 100)$pval, 0.05)
  expect_lt(test_multiplicity(Ysample, mult = c(2,2,2,1), 100)$pval, 0.05)
  expect_lt(test_multiplicity(Ysample, mult = c(4,1,1,1), 100)$pval, 0.05)
  expect_lt(test_multiplicity(Ysample, mult = c(3,1,3), 100)$pval, 0.05)
  
  # but power varies  
  expect_lt(test_multiplicity(Ysample, mult = c(3,2,2), 100)$pval, 0.3)
})

test_that("xiget() behaves properly", {
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
  n = 100
  mult <- c(3, 2)
  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  
  #if evecs supplied, should they be sorted so that the calculated approximate eigenvalues are descending?
  simxi <- function(n, mn, sigma, mult, idxs, evecs = NULL, sortevecs = FALSE){
    Ysample <- rsymm(n, mn, sigma)
    Ybar <- mmean(Ysample)
    if (is.null(evecs)){evals <- eigen_desc(Ybar)$values} 
    else {
      evals <- diag(t(evecs) %*% Ybar %*% evecs) #use the true eigenvectors
      if (sortevecs){
        ord <- order(evals, decreasing = TRUE)
        evals <- evals[ord]
      }
    }
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
  thecov <- xicovar(mult, idxs, eigen_desc(mn)$vectors, C0/n)
 
  set.seed(35468) 
  # semi-plugged in xi
  emcov_semi <- replicate(1000,
   simxi(n, mn = mn, sigma = C0, mult, idxs, eigen_desc(mn)$vectors)) |>
    t() |>
    cov()
  expect_equal(emcov_semi, thecov, tolerance = 0.01)
  
  # semi-plugged in xi with sorting: covariance changes a lot
  set.seed(35468) 
  emcov_sort <- replicate(1000,
    simxi(n, mn = mn, sigma = C0, mult, idxs, evecs = eigen_desc(mn)$vectors, sortevecs = TRUE)) |>
    t() |>
    cov()
  expect_equal(emcov_sort, thecov, tolerance = 0.01)
  
  set.seed(35468) 
  # fully empirical (plugged in) xi
  emcov <- replicate(1000,
    simxi(n, mn = mn, sigma = C0, mult, idxs)) |>
    t() |>
    cov()
  expect_equal(emcov, emcov_sort, tolerance = 0.01)
  expect_equal(emcov, thecov, tolerance = 0.01) # currently failing.
  # That is the covariance of fully empirical \hat\xi is similar to that of \xi with sorting
  # but vastly different to the covariance of \xi without sorting.
  # skip("Empirical covariance of xi using estimated evecs is very different to xi using population-mean evecs.")
  # expect_equal(emcov / thecov, matrix(1, 3, 3), tolerance = 0.5)
  #Empirical covariance here stable up to 100000 replicates, suggests the empirical xi (eigenvectors from the sample mean) has a very different distribution to the exact xi for n = 10, 100. Asymptotically, I'd expect it to be good though.
  #I wonder if changes in order of vectors are washing things out! But how to avoid that??
})

test_that("test p value resistant to fixed trace by normalisation", {
  # division by sum(vec) really changes the sizes of vec, so
  # the covariance of vec changes, so I'd expect the stat to be different
  # p values should be similar
  set.seed(134)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(10, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2))
  Ysample_n <- lapply(Ysample, normtrace)
  
  set.seed(34641)
  pval <- test_multiplicity(Ysample, mult = mult, 1000)$pval
  set.seed(34641)
  pval_n <- test_multiplicity(Ysample_n, mult = mult, 1000)$pval
  expect_equal(pval, 
               pval_n, tol = 1E-2)
})

test_that("fixed trace from projection preserved by standardisation and ignored by stat", {
  set.seed(134)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(1000, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2))
  Ysample_n <- lapply(Ysample, projtrace)
  
  # check that only first element of Hevals  is changed by trace fix
  evals <- eigen_desc(mmean(Ysample))$values
  evals_n <- eigen_desc(mmean(Ysample_n))$values
  expect_equal((helmert(sum(mult)) %*% evals)[-1],
    (helmert(sum(mult)) %*% evals_n)[-1])
  
  std <- standardise_multiplicity(Ysample, mult)
  std_n <- standardise_multiplicity(Ysample_n, mult)
  expect_true(all.equal(lapply(std, projtrace), std_n, check.attributes = FALSE))
  
  expect_equal(stat_multiplicity(Ysample, mult = mult), stat_multiplicity(Ysample_n, mult = mult))
})
