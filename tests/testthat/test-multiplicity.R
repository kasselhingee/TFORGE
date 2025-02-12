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
  mult <- c(3, 2, 1, 1)
  nullmean <- multiplicity_nullmean(av, mult = mult)
  # nullmean conserves trace:
  expect_equal(sum(diag(nullmean)), sum(diag(av)))
  
  # null mean has the same eigenvectors/space
  expect_equal(t(es$vectors) %*% nullmean %*% es$vectors, 
               diag(multiplicity_blk(es$values, mult = mult)))
  
  Ystdsample <- standardise_multiplicity(Ysample, mult = mult)
  
  # old way of calculating Ystdsample matches new way
  expect_equal(standardise_multiplicity_old(Ysample, mult = mult), Ystdsample)

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen_desc(newav)
  expect_equal(newes$values[1:3], rep(3, 3), tolerance = 0.1)
  expect_equal(newes$values[4:5], rep(2, 2), tolerance = 0.1)
  expect_equal(newes$values[6], es$values[6])
  expect_equal(newes$values[7], es$values[7])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 1, 1)), 0)
  expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 1, 1), refbasis = diag(1, 7)), 0)
  expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 1, 1), refbasis = runifortho(7)), 0)
  expect_error(expect_equal(stat_multiplicity(Ystdsample, mult = c(3, 2, 2)), 0))
  expect_error(expect_equal(stat_multiplicity(Ysample, mult = c(3, 2, 1, 1)), 0))
})

test_that("debugging stat with true evecs has correct null distribution", {
  set.seed(1332)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- replicate(1000, {
    Ysample <- rsymm_norm(200, diag(evals), sigma = diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    suppressWarnings(stat_multiplicity(Ysample, mult = mult, evecs = diag(sum(mult))))
  })
  
  # qqplot(vals, y = rchisq(1000, df = sum(mult-1)))
  res <- ks.test(vals, "pchisq", df = sum(mult-1))
  expect_gt(res$p.value, 0.1)
})

test_that("stat has correct null distribution", {
  set.seed(4)
  abasis <- runifortho(7)
  set.seed(1332)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- replicate(1000, {
    Ysample <- rsymm_norm(100, diag(evals), sigma = diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    c(statr = stat_multiplicity(Ysample, mult = mult),
      statc = stat_multiplicity(Ysample, mult = mult, refbasis = diag(1, 7)),
      statm = stat_multiplicity(Ysample, mult = mult, refbasis = "mincorr"),
      stata = stat_multiplicity(Ysample, mult = mult, refbasis = abasis))
  })
  
  # qqplot(vals["statr", ], y = rchisq(1000, df = sum(mult-1)))
  # qqplot(vals["statc", ], y = rchisq(1000, df = sum(mult-1)))
  # qqplot(vals["statm", ], y = rchisq(1000, df = sum(mult-1)))
  # qqplot(vals["stata", ], y = rchisq(1000, df = sum(mult-1)))
  expect_gt(ks.test(vals["statr", ], "pchisq", df = sum(mult-1))$p.value, 0.2)
  expect_gt(ks.test(vals["statc", ], "pchisq", df = sum(mult-1))$p.value, 0.2)
  expect_gt(ks.test(vals["statm", ], "pchisq", df = sum(mult-1))$p.value, 0.2)
  expect_gt(ks.test(vals["stata", ], "pchisq", df = sum(mult-1))$p.value, 0.15)
})

test_that("test has uniform distribution", {
  skip_on_cran() #test very slow
  set.seed(4)
  abasis <- runifortho(7)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  set.seed(5)#set.seed(1331)
  vals <- pbapply::pbreplicate(1000, { #1000 for more thorough
    Ysample <- rsymm_norm(100, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    # B = 100 for more thorough
    c(r = test_multiplicity(Ysample, mult = mult, B = 100)$pval,
      c = test_multiplicity(Ysample, mult = mult, B = 100, refbasis = diag(1, 7))$pval,
      m = test_multiplicity(Ysample, mult = mult, B = 100, refbasis = "mincorr")$pval,
      a = test_multiplicity(Ysample, mult = mult, B = 100, refbasis = abasis)$pval)
  })
  
  # qqplot(vals["r", ], y = runif(1000))
  # qqplot(vals["c", ], y = runif(1000))
  # qqplot(vals["m", ], y = runif(1000))
  # qqplot(vals["a", ], y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals["r", ], "punif"))$p.value, 0.05) #above 0.2 if above two thoroughness measures taken
  expect_gt(suppressWarnings(ks.test(vals["c", ], "punif"))$p.value, 0.05) 
  expect_gt(suppressWarnings(ks.test(vals["m", ], "punif"))$p.value, 0.05) 
  expect_gt(suppressWarnings(ks.test(vals["a", ], "punif"))$p.value, 0.05)
})

test_that("using sample evecs NOT does not give uniform p-values", {
  skip_on_cran() #test very slow
  set.seed(4)
  abasis <- runifortho(7)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  set.seed(5)#set.seed(1331)
  vals <- pbapply::pbreplicate(100, { #1000 for more thorough
    Ysample <- rsymm_norm(100, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    evecs <- eigen_desc(mmean(Ysample))$vectors
    test_multiplicity(Ysample, mult = mult, B = 100, refbasis = evecs)$pval
  })
  
  qqplot(vals, y = runif(1000))
  expect_lt(suppressWarnings(ks.test(vals, "punif"))$p.value, 0.001)
})

test_that("chisq: test has uniform distribution", {
  set.seed(5)
  abasis <- runifortho(7)
  set.seed(3) #set.seed(1331)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- replicate(100, { #1000 for more thorough
    Ysample <- rsymm_norm(100, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    test_multiplicity(Ysample, mult = mult, B = "chisq")$pval #B = 100 for more thorough
    c(r = test_multiplicity(Ysample, mult = mult, B = "chisq")$pval,
      c = test_multiplicity(Ysample, mult = mult, B = "chisq", refbasis = diag(1, 7))$pval,
      a = test_multiplicity(Ysample, mult = mult, B = "chisq", refbasis = abasis)$pval)
  })
  
  # qqplot(vals["r", ], y = runif(1000))
  # qqplot(vals["c", ], y = runif(1000))
  # qqplot(vals["a", ], y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals["r", ], "punif"))$p.value, 0.15)
  expect_gt(suppressWarnings(ks.test(vals["c", ], "punif"))$p.value, 0.15) 
  expect_gt(suppressWarnings(ks.test(vals["a", ], "punif"))$p.value, 0.15)
})

test_that("test rejects some incorrect hypotheses", {
  set.seed(13321)
  Ysample <- rsymm(100, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)),
                   sigma = 0.1 * diag(7 * (7 + 1)/2))
  set.seed(3654)
  res <- test_multiplicity(Ysample, mult = c(3,2,1,1), 100)
  expect_gt(res$pval, 0.1)

  # set.seed(3542)
  # pvals_2311 <- replicate(100, test_multiplicity(Ysample, mult = c(2,3,1,1), 100)$pval, cl = 3)
  set.seed(3542)
  expect_lt(test_multiplicity(Ysample, mult = c(2,3,1,1), 100)$pval, 0.05)
  set.seed(35423) 
  expect_lt(test_multiplicity(Ysample, mult = c(2,2,2,1), 100)$pval, 0.05)
  set.seed(35424) 
  expect_lt(test_multiplicity(Ysample, mult = c(4,1,1,1), 100)$pval, 0.05)
  set.seed(35425) 
  expect_lt(test_multiplicity(Ysample, mult = c(3,1,3), 100)$pval, 0.05)
  set.seed(35426) 
  expect_lt(test_multiplicity(Ysample, mult = c(3,3,1), 100)$pval, 0.05)
  
  set.seed(35427) 
  expect_lt(test_multiplicity(Ysample, mult = c(2,3,2), 100)$pval, 0.05)
  set.seed(3541) 
  expect_lt(test_multiplicity(Ysample, mult = c(3,2,2), 100)$pval, 0.05)
  # note that at B=100 there is still a lot a randomness in the output pvalue
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

test_that("xicovar() gives the same as sample covariance of xi", {
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
  C0_U <- runifortho(nvar)
  C0 <- C0_U %*% diag(1:nrow(C0_U)) %*% t(C0_U) #use this to simulate
  set.seed(348)
  mn_U <- runifortho(p)
  mn <- mn_U %*% diag(c(3, 3, 3, 2, 2)) %*% t(mn_U)
  
  # theoretical covariance
  thecov <- xicovar(mult, idxs, eigen_desc(mn)$vectors, C0/n)
 
  set.seed(35468) 
  # semi-empirical xi
  emcov_semi <- replicate(1E3,
   simxi(n, mn = mn, sigma = C0, mult, idxs, eigen_desc(mn)$vectors)) |>
    t() |>
    cov()
  expect_equal(emcov_semi, thecov, tolerance = 0.02)
  # note that xihat (eigenvectors given by mn) has a very different distribution because for each sample the eigenvectors for an eigenvalue can be any from the space,
  # not necessarily the exact ones supplied to thecov.
})

test_that("test p value resistant to fixed trace by normalisation", {
  # division by sum(vec) really changes the sizes of vec, so
  # the covariance of vec changes, so I'd expect the stat to be different
  # p values should be similar
  set.seed(134)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(10, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2))
  Ysample_n <- normalise_trace(Ysample)
  
  set.seed(34641)
  pval <- test_multiplicity(Ysample, mult = mult, 1000, refbasis = diag(1, 7))$pval
  set.seed(34641)
  pval_n <- test_multiplicity(Ysample_n, mult = mult, 1000, refbasis = diag(1, 7))$pval
  expect_equal(pval, 
               pval_n, tol = 1E-1)
})

test_that("runifortho produces orthogonal matrices", {
  set.seed(354)
  rot <- runifortho(3)
  expect_equal(rot %*% t(rot), diag(3))
  set.seed(354796)
  rot <- runifortho(6)
  expect_equal(rot %*% t(rot), diag(6))
})

test_that("runifortho density is invariant to rotations", {
  set.seed(354)
  rot <- runifortho(3)
  set.seed(3514)
  mats1 <- replicate(1E3, runifortho(3))
  mats2 <- replicate(1E3, runifortho(3)) #for full independencs
  # mats <- replicate(1E3, rstiefel::rustiefel(3,3))
  
  #get number of mats that have last column within certain amount of e1=c(1,0,0)
  orig_in <- mats1[1,3,] < pi/16
  origest <- list(
    mean = mean(orig_in),
    se = sd(orig_in)/sqrt(length(orig_in))
  )
  
  rotmats <- simplify2array(apply(mats2, 3, FUN = '%*%', rot, simplify = FALSE))
  rot_in <- rotmats[1,3,]  < pi/16
  rotest <- list(
    mean = mean(rot_in),
    se = sd(rot_in)/sqrt(length(rot_in))
  )
  expect_lt(abs(rotest$mean - origest$mean), 2 * sqrt(rotest$se^2 + origest$se^2))
})
