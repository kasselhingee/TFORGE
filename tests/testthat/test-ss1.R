test_that("stat_ss1() on single sample from NULL is consistent with chisq", {
  set.seed(13)
  vals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones
        evecs <- eigen(m)$vectors
        evals <- eigen(m)$values
        evals <- evals/sqrt(sum(evals^2))
        out <- evecs %*% diag(evals) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        return(out)
    })
    stat_ss1(Y, evals = c(3,2,1)/sqrt(sum(c(3,2,1)^2)))
    })
  
  # qqplot(vals, y = rchisq(1000, df = 2))
  res <- ks.test(vals, "pchisq", df = 2)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1() on multiple NULL samples is consistent with chisq", {
  set.seed(13)
  vals <- replicate(100, {
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones
      evecs <- eigen(m)$vectors
      evals <- eigen(m)$values
      evals <- evals/sqrt(sum(evals^2))
      out <- evecs %*% diag(evals) %*% t(evecs)
      out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
      return(out)
    })
  }, simplify = FALSE)
  stat_ss1(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*2))
  res <- ks.test(vals, "pchisq", df = (5-1)*2)
  expect_gt(res$p.value, 0.2)
})

test_that("test_ss1() Omega from sim evals in sst matches", {
  set.seed(1221)
  evals <- mvtnorm::rmvnorm(50, mean = c(1, 0, 0))
  # evals <- evals / sqrt(rowSums(evals^2))
  # pick a fixed set of vectors randomly
  Q <- eigen(rsymm_norm(1, mean = diag(c(3,2,1)))[[1]])$vectors
  Y <- apply(evals, 1, function(v){t(Q) %*% diag(v) %*% Q}, simplify = FALSE)
  expect_equal(cov_evals(Y), cov(evals))
  
  #expect that projection of covariance onto the tangent at c(1,0,0) will 
  #keep cov(evals) the same for the second two directions, and zero for the first dimension
  Delta <- amaral2007Lemma1(c(1, 0, 0)/sqrt(sum(c(1, 0, 0)^2)))
  Delta %*% cov(evals) %*% t(Delta)
  stat_ss1(Y, evals = c(1,0,0))
  res <- test_ss1(Y, c(1,0,0), 100, maxit = 25)
  res$pval
  
})

test_that("test_ss1() pval on NULL Normal evals sst is uniform", {
  set.seed(1333)
  pvals <- replicate(100, {
    evals <- mvtnorm::rmvnorm(50, mean = c(3, 2, 1))
    # evals <- evals / sqrt(rowSums(evals^2))
    # pick a fixed set of vectors randomly
    Q <- eigen(rsymm_norm(1, mean = diag(c(3,2,1)))[[1]])$vectors
    Y <- apply(evals, 1, function(v){t(Q) %*% diag(v) %*% Q}, simplify = FALSE)
    Y <- lapply(Y, function(m){m[lower.tri(m)] <- m[upper.tri(m)]; m}) #remove machine error
    res <- test_ss1(Y, c(3,2,1), 100, maxit = 100)
    res$pval
  })
  qqplot(pvals, y = runif(100))
  expect_gt(ks.test(pvals, "punif")$p.value, 0.05)
})

test_that("test_ss1() pval on NULL Normal sst is uniform", {
  set.seed(1333)
  pvals <- replicate(100, {
    Y <- rsymm_norm(300, diag(c(3,2,1)), sigma = diag(1, 6)) #at 300 samples it works :)
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones
      evecs <- eigen(m)$vectors
      evals <- eigen(m)$values
      evals <- evals/sqrt(sum(evals^2))
      out <- evecs %*% diag(evals) %*% t(evecs)
      out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
      return(out)
    })
    res <- test_ss1(Y, c(3,2,1), 100, maxit = 25)
    res$pval
  })
  qqplot(pvals, y = runif(100))
  expect_gt(ks.test(pvals, "punif")$p.value, 0.05)
})

test_that("test_ss1() reject for single sample with wrong eval", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(3,2,1)), sigma = diag(0.7, 6))
  Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones
    evecs <- eigen(m)$vectors
    evals <- eigen(m)$values
    evals <- evals/sqrt(sum(evals^2))
    out <- evecs %*% diag(evals) %*% t(evecs)
    out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
    return(out)
  })

  res <- test_ss1(Y, c(1,1,1)/3, 100, maxit = 100)
  expect_lt(res$pval, 0.05)
})


test_that("amaral2007Lemma1() produces correct result for a unit vector", {
  m <- runif(5, -1, 1)
  m <- m/sqrt(sum(m^2))
  
  A <- amaral2007Lemma1(m)
  expect_equal(drop(A %*% m), rep(0, 4))
  expect_equal(A %*% Conj(t(A)), diag(1, 4))
})
