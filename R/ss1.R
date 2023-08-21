#' @title Methods for testing eigenvalues with sum of squares = 1
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
stat_ss1 <- function(x, evals = NULL, evecs = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  if (is.null(evals) && inherits(x, "sst")){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  if (!is.null(evecs) && (length(x) > 1)){warning("evecs supplied for multisample situation supplied. This is unusual.")}
  
  # means and eigenspaces used in multiple parts, so calculate first here:
  mns <- lapply(x, mmean)
  ess <- lapply(mns, eigen)
  evalsav <- lapply(ess, "[[", "values")
  
  # first get each Omega2j:
  covars_unconstrained <- mapply(cov_evals, ms = x, evecs = lapply(ess, "[[", "vectors"), av = mns, SIMPLIFY = FALSE)
  Deltas <- lapply(evalsav, function(d2){amaral2007Lemma1(d2/sqrt(sum(d2^2)))})
  Omega2s <- mapply(function(d2, covar, Delta){
    Delta %*% covar %*% t(Delta) / sum(d2^2)
  }, d2 = evalsav, covar = covars_unconstrained, Delta = Deltas, SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
    #estimate according to (36)
    mats <- mapply(function(Delta, Omega){t(Delta) %*% solve_NAonerror(Omega, NAonerror) %*% Delta},
                   Delta = Deltas,
                   Omega = Omega2s, SIMPLIFY = FALSE)
    mat <- purrr::reduce(mats, `+`)
    d0 <- eigen(mat)$vectors[, ncol(mat)]
  } else {
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
  }
  
  # now the statistic (32) for each sample:
  persamplestat <- mapply(function(d3, Delta, Omega, n){
    n * t(d3/sqrt(sum(d3^2)) - d0) %*% t(Delta) %*% solve_NAonerror(Omega, NAonerror) %*% Delta %*% (d3/sqrt(sum(d3^2)) - d0)
  },
  d3 = evalsav, #not yet normalised as in (32)
  Delta = Deltas,
  Omega = Omega2s,
  n = lapply(x, length),
  SIMPLIFY = FALSE
  )
  stat <- drop(purrr::reduce(persamplestat, `+`))
  attr(stat, "null_evals") <- drop(d0)
  return(stat)
}

# the construction of matrix A in Lemma1 Amaral, G. J. A., Dryden, I. L., & Wood, A. T. A. (2007). Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. Journal of the American Statistical Association, 102(478), 695–707. http://www.jstor.org/stable/27639898
 
amaral2007Lemma1 <- function(m){
  d <- length(m)
  finalelement <- m[d]
  # the following is a rearrangement to make result less computationally sensitive to size of final element.
  A1 <- diag(1, d-1) - m[-d] %*% Conj(t(m[-d]))/(1+abs(finalelement))
  if (finalelement < 0){
    A1 <- -A1
  }
  A <- cbind(A1, -m[-d])
  return(A)
}

# wrapper around solve that returns a matrix of NA if couldn't solve
solve_NAonerror <- function(A, NAonerror){
  erroraction <- function(e){
    if (!NAonerror){stop(e)}
    else {
      out <- NA*A
      attr(out, "message") <- e$message
      return(out)
    }
  }
  out <- tryCatch(solve(A), error = erroraction)
  out
}

#' @describeIn stat_ss1 Bootstrap test.
#' @param maxit The maximum number of iterations to use in finding the weights. Passed to `[emplik::el.test()]`.
#' @details The test did not perform well when the dispersion was very high (e.g. Normal entries with mean diagonal c(3,2,1) and variance of 1) - more studying needed.
#' @export
test_ss1 <- function(mss, evals = NULL, B, maxit = 25){
  mss <- as.mstorsst(mss)
  if (inherits(mss, "sst")){mss <- as.mstorsst(list(mss))}
  t0 <- stat_ss1(mss, evals = evals, NAonerror = FALSE)
  evals <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to evals)
  # also compute the bounds on possible cj in equation (37). See Eq37_cj_bound.pdf
  nullmeans <- lapply(mss, function(ms){
    av <- mmean(ms)
    evecs <- eigen(av)$vectors
    nullmean <- evecs %*% diag(evals) %*% t(evecs)
    # bounds for cj
    diags <- lapply(ms, function(m){diag(t(evecs) %*% m %*% evecs)})
    diags <- do.call(cbind, diags)
    ranges <- t(apply(diags, 1, range))
    colnames(ranges) <- c("min", "max")
    # incorporate the proposed eigenvalues:
    ranges <- ranges/evals
    ranges <- t(apply(ranges, 1, sort))
    # take largest min and smallest max as range
    crange <- c(min = max(ranges[, 1]), max = min(ranges[, 2]))
    
    # return
    attr(nullmean, "c_range") <- crange
    return(nullmean)
  })
  
  # compute corresponding weights that lead to emp.lik.
  # note that the profile likelihood function (result of el.test) has convex superlevel sets according to Theorem 3.2 (Owen 2001).
  # So there is unique minimum value to the problem where the mean lies on a line.
  wts <- mapply(function(ms, nullmean){
    if (attr(nullmean, "c_range")[["min"]] > attr(nullmean, "c_range")[["max"]]){
      return(rep(0, length(ms))) #no pluasible values of c - avoids error triggered in optimise
    }
    msarr <- do.call(rbind, lapply(ms, vech))
    bestmult <- optimise(f = function(x){#optim warns that Nelder-Mead unreliable on 1 dimension so using Brent here instead
            elres <- emplik::el.test(msarr, x*vech(nullmean), maxit = maxit)
            return(elres[["-2LLR"]])
          },
          lower = attr(nullmean, "c_range")[["min"]], 
          upper = attr(nullmean, "c_range")[["max"]]) 
    elres <- emplik::el.test(msarr, bestmult$minimum*vech(nullmean), maxit = maxit)
    if (elres$nits == maxit){warning("Reached maximum iterations in el.test() at best null mean.")}
    elres$wts
  }, ms = mss, nullmean = nullmeans, SIMPLIFY = FALSE)
  
  #check the weights
  wtsums_discrepacies <- vapply(wts, function(x){abs(length(x) - sum(x))}, FUN.VALUE = 0.1)
  if (any(wtsums_discrepacies > 1E-2)){
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is close or outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = Inf,
      nullt = NA,
      B = NA
    ))
  }
  
  res <- bootresampling(mss, wts, 
                        stat = stat_ss1,
                        B = B)
  return(res)
}