test_that("stat_ss1fixedtrace() on single sample from NULL is consistent with chisq", {
  set.seed(1311)
  vals <- pbreplicate(2000, {
    Y <- rsymm_norm(30, diag(c(1, 0, -1)))
    Y <- lapply(Y, function(m) {diag(m) <- diag(m) - drop(diag(m) %*% rep(1/sqrt(3), 3)) * rep(1/sqrt(3), 3); return(m)})
    Y <- lapply(Y, function(m) { #replace eigenvalues with normalised ones
        evecs <- eigen(m)$vectors
        evals <- eigen(m)$values
        evals <- evals/sqrt(sum(evals^2))
        out <- evecs %*% diag(evals) %*% t(evecs)
        out[lower.tri(out)] <- out[upper.tri(out)] #to remove machine differences
        return(out)
    })
    # hasfixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps))
    # hasss1(Y)
    stat_ss1fixedtrace(Y, evals = c(1, 0, -1))
    }, cl = 2)
  
  qqplot(vals, y = rchisq(1E6, df = 1))
  res <- ks.test(vals, "pchisq", df = 1)
  expect_gt(res$p.value, 0.2)
})