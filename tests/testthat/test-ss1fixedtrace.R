test_that("stat_ss1fixedtrace() on single sample from NULL with fixed evecs is not inconsistent with chisq", {
  skip_if_fast_check() #non-core functionality
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
    evecs <- eigen_desc(invvech(rsymm_norm(1, mean = diag(1, 3))[1, ]))$vectors
    evals <- revals(50, m = c(1/sqrt(2), 0, -1/sqrt(2)))
    Y <- apply(evals, 1, function(v){
      out <- evecs %*% diag(v) %*% t(evecs)
      out <- makeSymmetric(out) #to remove machine differences
      out
      }, simplify = FALSE)
    Y <- as_fsm(Y)
    stat_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)))
  })
  # qqplot(vals, y = rchisq(1E6, df = 1))
  res <- ks.test(vals, "pchisq", df = 1)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1fixedtrace() on single sample from NULL is consistent with chisq", {
  set.seed(1311)
  vals <- replicate(ifelse(fast_check_on(), 100, 1000), {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    # has_fixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps))
    # has_ss1(Y)
    stat_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)))
    })
  
  # qqplot(vals, y = rchisq(1E6, df = 1))
  res <- ks.test(vals, "pchisq", df = 1)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_ss1fixedtrace() on TFORGE_kfsm from NULL is not inconsistent with chisq on n=300", {
  vals <- vapply(300 + (1:ifelse(fast_check_on(), 100, 1000)), function(seed){
    set.seed(seed)
    Ysamples <- lapply(c(1000, 300), function(n){
    Y <- rsymm_norm(n, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    Y})
    stat_ss1fixedtrace(Ysamples)
    }, FUN.VALUE = 1.32)
  
  # qqplot(vals, y = rchisq(1E6, df = (2-1) * 1))
  res <- ks.test(vals, "pchisq", df = (2-1) * 1)
  expect_gt(res$p.value, 0.2)
})

test_that("test_ss1fixedtrace() uniform pval on NULL TFORGE_fsm", {
  set.seed(1333)
  pvals <- replicate(ifelse(fast_check_on(), 20, 100), {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    stopifnot(has_fixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps)))
    stopifnot(has_ss1(Y))
    res <- test_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)), B = ifelse(fast_check_on(), 20, 100), maxit = 1000)
    res$pval
    })
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("chisq: test_ss1fixedtrace() uniform pval on NULL TFORGE_fsm", {
  set.seed(1333)
  pvals <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
    Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    stopifnot(has_fixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps)))
    stopifnot(has_ss1(Y))
    res <- test_ss1fixedtrace(Y, evals = c(1/sqrt(2), 0, -1/sqrt(2)), B = "chisq", maxit = 1000)
    res$pval
  })
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test_ss1fixedtrace() uniform pval on NULL TFORGE_kfsm", {
  pvals <- vapply(13131 + (1:ifelse(fast_check_on(), 20, 100)), function(seed){
    set.seed(seed)
    Ysamples <- replicate(2, {
      Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
      Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
      Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    }, simplify = FALSE)
    res <- test_ss1fixedtrace(Ysamples, B = ifelse(fast_check_on(), 20, 100), maxit = 1000)
    res$pval
  }, FUN.VALUE = 1.3)
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.2)
})

test_that("chisq: test_ss1fixedtrace() uniform pval on NULL TFORGE_kfsm", {
  pvals <- vapply(13131 + (1:100), function(seed){
    set.seed(seed)
    Ysamples <- replicate(2, {
      Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))))
      Y <- project_trace(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
      Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    }, simplify = FALSE)
    res <- test_ss1fixedtrace(Ysamples, B = "chisq", maxit = 1000)
    res$pval
  }, FUN.VALUE = 1.3)
  # qqplot(pvals, y = runif(100))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.2)
})
