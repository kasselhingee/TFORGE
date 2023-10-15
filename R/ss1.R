#' @title Methods for testing eigenvalues with sum of squares = 1
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
stat_ss1 <- function(x, evals = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  
  # means and eigenspaces used in multiple parts, so calculate first here:
  mns <- lapply(x, mmean)
  ess <- lapply(mns, eigen_desc)
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
    d0 <- eigen_desc(mat)$vectors[, ncol(mat)]
    # make d0 DOT evalsav have as much positive sign as possible
    dotprds <- vapply(evalsav, function(v){v %*% d0}, FUN.VALUE = 0.1)
    avsign <- mean(sign(dotprds))
    if (avsign < 0){
      d0 <- -1 * d0
    }
    if (!all(order(d0, decreasing = TRUE) == 1:length(d0))){warning("Estimated common eigenvalues are not in decreasing order.")}
  } else {
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
  }
  
  # now the statistic (32) for each sample:
  persamplestat <- mapply(function(d2, Delta, Omega, n){
    n * t(d2/sqrt(sum(d2^2)) - d0) %*% t(Delta) %*% solve_NAonerror(Omega, NAonerror) %*% Delta %*% (d2/sqrt(sum(d2^2)) - d0)
  },
  d2 = evalsav, #not yet normalised as in (32)
  Delta = Deltas,
  Omega = Omega2s,
  n = lapply(x, nrow),
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
test_ss1 <- function(mss, evals = NULL, B, maxit = 25, sc = TRUE){
  mss <- as.mstorsst(mss)
  stopifnot(hasss1(mss))
  if (inherits(mss, "sst")){mss <- as.mstorsst(list(mss))}
  if (is.null(evals) && (length(mss) == 1)){stop("evals must be supplied for a meaningful test since mss is a single sample")}
  if (!is.null(evals) && (length(mss) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  
  t0 <- stat_ss1(mss, evals = evals, NAonerror = FALSE)
  d0 <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to d0)
  # also compute the bounds on possible cj in equation (37). See Eq37_cj_bound.pdf
  nullmeans <- lapply(mss, elnullmean, d0 = d0, getcbound = TRUE)
  
  # compute corresponding weights that lead to emp.lik.

  wts <- mapply(opt_el.test, ms = mss, mu = nullmeans, maxit = maxit, sc = sc, SIMPLIFY = FALSE)
  
  #check the weights
  wtsums_discrepacies <- vapply(wts, function(x){abs(length(x) - sum(x))}, FUN.VALUE = 0.1)
  if (any(wtsums_discrepacies > 1E-2)){
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is close or outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    ))
  }
  
  res <- bootresampling(mss, wts, 
                        stat = stat_ss1,
                        B = B,
                        evals = evals)
  return(res)
}

#' @title Test whether the supplied sample(s) have ss1
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @export
#' @return `TRUE` or `FALSE`
hasss1 <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as.mstorsst(x)
  if (inherits(x, "mst")){x <- do.call(rbind, x)}
  ss <- apply(x, 1, function(v){
    m <- invvech(v)
    sum(eigen_desc(m)$values^2)
  })
  isTRUE(all.equal(ss, rep(1, length(ss)), tolerance = tolerance))
}


# compute means that satisfy the NULL hypothesis (eigenvalues equal to d0) for use in empirical likelihood
# also compute the bounds on possible cj in equation (37). See Eq37_cj_bound.pdf
# @param a single sample of tensors
# @param d0 the NULL set of eigenvalues
# @param av The average of `ms`, if omitted then computed from `ms` directly (include to save computation time)
# @param evecs The eigenvectors of the average of `ms`. If omitted then computed from the average of `ms` directly (include to save computation time).
elnullmean <- function(ms, d0, av = NULL, evecs = NULL, getcbound = FALSE){
  if (is.null(av)){av <- mmean(ms)}
  if (is.null(evecs)){evecs <- eigen_desc(av)$vectors}
  nullmean <- evecs %*% diag(d0) %*% t(evecs)
  if (getcbound){# bounds for cj
    diags <- apply(ms, 1, function(vec){m <- invvech(vec); diag(t(evecs) %*% m %*% evecs)}, simplify = FALSE)
    diags <- do.call(cbind, diags)
    ranges <- t(apply(diags, 1, range))
    colnames(ranges) <- c("min", "max")
    # incorporate the proposed eigenvalues:
    ranges <- ranges/d0
    ranges <- t(apply(ranges, 1, sort))
    # take largest min and smallest max as range
    crange <- c(min = max(ranges[, 1]), max = min(ranges[, 2]))
    
    attr(nullmean, "c_range") <- crange
  }
  return(nullmean)
}

#function finds the best c and weights such that weighted average of data is c.mu and 
#empirical likelihood maximised
# note that the profile likelihood function (result of el.test) has convex superlevel sets according to Theorem 3.2 (Owen 2001).
# So there is unique minimum value to the problem where the mean lies on a line.
# @param ms A single sample of symmetric tensors
# @param mu Proposed mean up-to-constant c. It is assumed to have an attribute "c_range" range gives a range of values of c, passed to `optimize()`
# @param sc If TRUE use Owen's self-concordant emplik function
opt_el.test <- function(ms, mu, maxit = 25, sc = FALSE){
  class(ms) <- "matrix"
  if (attr(mu, "c_range")[["min"]] > attr(mu, "c_range")[["max"]]){
    return(rep(0, nrow(ms))) #no pluasible values of c - avoids error triggered in optimise
  }
  if (sc) {
    bestmult <- optimise(f = function(x){#optim warns that Nelder-Mead unreliable on 1 dimension so using Brent here instead
      scelres <- emplik(ms, x*vech(mu), itermax = maxit)
      return(-scelres$logelr)
    },
    lower = attr(mu, "c_range")[["min"]], 
    upper = attr(mu, "c_range")[["max"]]) 
    scelres <- emplik(ms, bestmult$minimum*vech(mu), itermax = maxit)
    if (!isTRUE(scelres$converged)){warning("emplik() did not converge")}
    wts <- as.vector(scelres$wts) * nrow(ms) #the multiple here is to match the weights put out by emplik::el.test()
    # check result with el.test
    elres <- emplik::el.test(ms, 
                             bestmult$minimum*vech(mu), 
                             lam = as.vector(scelres$lam),
                             maxit = maxit)
    if (!isTRUE(all.equal(elres$wts, wts))){
      warning("Weights from emplik() differ from check by emplik::el.test()")
    }
  } else {
    bestmult <- optimise(f = function(x){#optim warns that Nelder-Mead unreliable on 1 dimension so using Brent here instead
            elres <- emplik::el.test(ms, x*vech(mu), maxit = maxit)
            return(elres[["-2LLR"]])
          },
          lower = attr(mu, "c_range")[["min"]], 
          upper = attr(mu, "c_range")[["max"]]) 
    elres <- emplik::el.test(ms, bestmult$minimum*vech(mu), maxit = maxit)
    if (elres$nits == maxit){warning("el.test() reached maximum iterations of ", maxit, " at best null mean.")}
    wts <- elres$wts
  }
  wts
}


#' Eigenvalues divided to have squared sum 1 
#' Uses the fact that A*A squares the eigenvalues of A, and the trace of a matrix is the sum of the eigenvalues.
#' m A symmetric matric
normL2evals <- function(m){
  ssq <- sum(diag(m %*% m))
  newm <- m / sqrt(ssq)
  newm <- makeSymmetric(newm) #remove computational inaccuracies
  return(newm)
}
# ms is an sst
normL2evals_sst <- function(ms){
  as.sst(apply(ms, 1, function(v){normL2evals(invvech(v))}, simplify = FALSE))
}
