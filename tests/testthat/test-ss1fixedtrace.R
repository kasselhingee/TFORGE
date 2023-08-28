test_that("stat_ss1fixedtrace() on single sample from NULL with fixed evecs is not inconsistent with chisq", {
  set.seed(1353)
  #method for simulating eigenvalues
  revals <- function(n, m = c(1/sqrt(2), 0, -1/sqrt(2))){
    H <- helmertsub(3)
    projevals <- mvtnorm::rmvnorm(n, mean = H %*% m)
    projevals <- projevals/sqrt(rowSums(projevals^2))
    evals <- projevals %*% H
    evals
  }
  vals <- replicate(1000, {
    # choose a fixed set of eigenvectors
    evecs <- eigen(rsymm_norm(1, mean = diag(1, 3))[[1]])$vectors
    evals <- revals(50, m = c(1/sqrt(2), 0, -1/sqrt(2)))
    Y <- apply(evals, 1, function(v){
      out <- evecs %*% diag(v) %*% t(evecs)
      out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
      out
      }, simplify = FALSE)
    stat_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)))
  })
  # qqplot(vals, y = rchisq(1E6, df = 1))
  res <- ks.test(vals, "pchisq", df = 1)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1fixedtrace() on single sample from NULL is consistent with chisq", {
  set.seed(1311)
  vals <- replicate(1000, {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)}) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
        evecs <- eigen(m)$vectors
        evals <- eigen(m)$values
        evals <- evals/sqrt(sum(evals^2))
        out <- evecs %*% diag(evals) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        return(out)
    })
    # hasfixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps))
    # hasss1(Y)
    stat_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)))
    })
  
  # qqplot(vals, y = rchisq(1E6, df = 1))
  res <- ks.test(vals, "pchisq", df = 1)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1fixedtrace() on mst from NULL with fixed evecs per sample is not inconsistent with chisq", {
  set.seed(1353351)
  #method for simulating eigenvalues
  revals <- function(n, m = c(1/sqrt(2), 0, -1/sqrt(2))){
    H <- helmertsub(3)
    projevals <- mvtnorm::rmvnorm(n, mean = H %*% m)
    projevals <- projevals/sqrt(rowSums(projevals^2))
    evals <- projevals %*% H
    evals
  }
  vals <- replicate(1000, {
    Ysamples <- replicate(2, {
      # choose a fixed set of eigenvectors
      evecs <- eigen(rsymm_norm(1, mean = diag(1, 3))[[1]])$vectors
      evals <- revals(50, m = c(1/sqrt(2), 0, -1/sqrt(2)))
      Y <- apply(evals, 1, function(v){
        out <- evecs %*% diag(v) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        out
      }, simplify = FALSE)
      Y
    }, simplify = FALSE)
    stat_ss1fixedtrace(Ysamples)
  })
  # qqplot(vals, y = rchisq(1E6, df = (2-1) * 1))
  res <- ks.test(vals, "pchisq", df = (2-1) * 1)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1fixedtrace() on mst from NULL is not inconsistent with chisq on n=200", {
  vals <- vapply(300 + (1:1000), function(seed){
    set.seed(seed)
    Ysamples <- replicate(2, {
    Y <- rsymm_norm(300, diag(c(1/sqrt(2), 0, -1/sqrt(2)))) #samples of 300 are big enough, but 200 are not
    Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)}) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
        evecs <- eigen(m)$vectors
        evals <- eigen(m)$values
        evals <- evals/sqrt(sum(evals^2))
        out <- evecs %*% diag(evals) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        return(out)
    })
    Y}, simplify = FALSE)
    stat_ss1fixedtrace(Ysamples)
    }, FUN.VALUE = 1.32)
  
  # qqplot(vals, y = rchisq(1E6, df = (2-1) * 1))
  res <- ks.test(vals, "pchisq", df = (2-1) * 1)
  expect_gt(res$p.value, 0.2)
})

test_that("test_ss1fixedtrace() uniform pval on NULL sst", {
  set.seed(1333)
  pvals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)}) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
        evecs <- eigen(m)$vectors
        evals <- eigen(m)$values
        evals <- evals/sqrt(sum(evals^2))
        out <- evecs %*% diag(evals) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        return(out)
    })
    stopifnot(hasfixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps)))
    stopifnot(hasss1(Y))
    res <- suppressWarnings(test_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)), B = 100, maxit = 1000))
    res$pval
    })
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test_ss1fixedtrace() uniform pval on NULL mst", {
  pvals <- vapply(13131 + (1:100), function(seed){
    set.seed(seed)
    Ysamples <- replicate(2, {
      Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
      Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)}) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
      Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
          evecs <- eigen(m)$vectors
          evals <- eigen(m)$values
          evals <- evals/sqrt(sum(evals^2))
          out <- evecs %*% diag(evals) %*% t(evecs)
          out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
          return(out)
      })
    }, simplify = FALSE)
    res <- suppressWarnings(test_ss1fixedtrace(Ysamples, B = 100, maxit = 1000))
    res$pval
  }, FUN.VALUE = 1.3)
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test_ss1fixedtrace() returns on bad el weights for mst", {
  set.seed(13133)
  Ysamples <- replicate(2, {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)}) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
      evecs <- eigen(m)$vectors
      evals <- eigen(m)$values
      evals <- evals/sqrt(sum(evals^2))
      out <- evecs %*% diag(evals) %*% t(evecs)
      out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
      return(out)
    })
  }, simplify = FALSE)
  suppressWarnings(res <- test_ss1fixedtrace(Ysamples, B = 100, maxit = 1000))
  expect_gt(length(res$nullt), 1)
  expect_lt(min(vapply(res$stdx, sum, FUN.VALUE = 1.3)), 49)
})

