#' @title Test for eigenvalues when sum of squared eigenvalues is 1
#' @description 
#' For a single sample of symmetric matrices where sum of squared eigenvalues  = 1, test eigenvalues of the population mean.
#' For multiple samples of symmetric matrices where sum of squared eigenvalues  = 1, test for equality of the eigenvalues of the population means.
#' The test statistic is calculated by `stat_ss1()`.
#' @details
#' The sum of squared eigenvalues constraint forces the set of eigenvalues to lie on a sphere (or circle).
#' The test statistic accounts for this constraint by projecting eigenvalues onto a plane perpendicular to the direction of the sample average's eigenvalues.
#' Bootstrap resampling is from an empirical distribution that satisfies the null hypothesis; for this test we use empirical likelihood \insertCite{owen:2013}{TFORGE} to find probability mass weights for each matrix in the original sample.
#'
#' Eigenvalues must be distinct.
#' @inheritParams test_unconstrained
#' @inheritParams test_fixedtrace
#' @inheritSection test_fixedtrace Hypotheses
#' @references \insertAllCited{}
#' @export
test_ss1 <- function(x, evals = NULL, B = 1000, maxit = 25){
  x <- as_flat(x)
  stopifnot(has_ss1(x))
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  
  if (has_fixedtrace(x)){warning("All tensors the same trace. Consider using test_ss1fixedtrace().")}
  
  t0 <- stat_ss1(x, evals = evals)
  d0 <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to d0)
  # also compute the bounds on possible cj in equation (37). See Eq37_cj_bound.pdf
  nullmeans <- lapply(x, elnullmean, d0 = d0, getcbound = TRUE)
  
  # compute corresponding weights that lead to emp.lik.

  wts <- mapply(opt_el.test, ms = x, mu = nullmeans, maxit = maxit, SIMPLIFY = FALSE)
  
  #check the weights
  if (!wtsokay(wts)){
    out <- list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    )
    class(out) <- c("TFORGE", class(out))
    return(out)
  }
  
  res <- bootresampling(x, wts, 
                        stat = stat_ss1,
                        B = B,
                        evals = evals)
  return(res)
}

stat_ss1 <- function(x, evals = NULL){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  
  # means and eigenspaces used in multiple parts, so calculate first here:
  mns <- lapply(x, mmean)
  ess <- lapply(mns, eigen_desc)
  evalsav <- lapply(ess, "[[", "values")
  
  # first get each Omega2j:
  covars_unconstrained <- mapply(cov_evals_est, ms = x, evecs = lapply(ess, "[[", "vectors"), av = mns, SIMPLIFY = FALSE)
  Deltas <- lapply(evalsav, function(d2){amaral2007Lemma1(d2/sqrt(sum(d2^2)))})
  Omega2s <- mapply(function(d2, covar, Delta){
    Delta %*% covar %*% t(Delta) / sum(d2^2)
  }, d2 = evalsav, covar = covars_unconstrained, Delta = Deltas, SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
    #estimate according to (35)
    mats <- mapply(function(n, Delta, Omega){n * t(Delta) %*% solve_error(Omega) %*% Delta},
                   n = vapply(x, nrow, FUN.VALUE = 1),
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
    if (!all(order(d0, decreasing = TRUE) == 1:length(d0))){
      d0 <- descendingordererror(d0)
    }
  } else {
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
  }
  
  # now the statistic (32) for each sample:
  persamplestat <- mapply(function(d2, Delta, Omega, n){
    n * t(d2/sqrt(sum(d2^2)) - d0) %*% t(Delta) %*% solve_error(Omega) %*% Delta %*% (d2/sqrt(sum(d2^2)) - d0)
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


#' @title Check whether the supplied sample(s) have equal sum of squared eigenvalues
#' @description Compares whether the sum of the squared eigenvalues of the supplied matrices match each other using the property that the sum of the squared eigenvalues of `Y` equals the trace of `Y %*% Y`.
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @export
#' @return `TRUE` or `FALSE`
has_ss1 <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_kfsm")){x <- do.call(rbind, x)}
  ss <- apply(x, 1, function(v){
    m <- invvech(v)
    sum(diag(m%*%m))
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
# note that the profile likelihood function (result of emplik()) has convex superlevel sets according to Theorem 3.2 (Owen 2001).
# So there is unique minimum value to the problem where the mean lies on a line.
# @returns If the method doesn't converge then the negative of the weights is returned
# @param ms A single sample of symmetric tensors
# @param mu Proposed mean up-to-constant c. It is assumed to have an attribute "c_range" range gives a range of values of c, passed to `optimize()`
opt_el.test <- function(ms, mu, maxit = 25){
  if (attr(mu, "c_range")[["min"]] > attr(mu, "c_range")[["max"]]){
    return(rep(0, nrow(ms))) #no pluasible values of c - avoids error triggered in optimise
  }
  bestmult <- optimise(f = function(x){#optim warns that Nelder-Mead unreliable on 1 dimension so using Brent here instead
    scelres <- emplik(ms, x*vech(mu), itermax = maxit)
    return(-scelres$logelr)
  },
  lower = attr(mu, "c_range")[["min"]], 
  upper = attr(mu, "c_range")[["max"]]) 
  scelres <- emplik(ms, bestmult$minimum*vech(mu), itermax = maxit)
  wts <- as.vector(scelres$wts) * nrow(ms) #the multiple here is to match the weights put out by emplik::el.test()

  # check and warn about results
  if (!isTRUE(scelres$converged)){
    wts <- -1 * wts
  }
  wts
}

# A function for checking that a list of weights is okay
wtsokay <- function(wts){
  if (any(vapply(wts, sum, FUN.VALUE = 1.3) < 0)){ #negative weights
    warning("The proposed null mean appears to be outside the convex hull for at least one sample (emplik() did not converge)" , call. = FALSE)
    return(FALSE)
  }
  if (any(vapply(wts, sum, FUN.VALUE = 1.3) < 0.5)){# weights sum to zero
    warning("The proposed null mean appears to be outside the convex hull for at least one sample (empirical likelihood weights are zero)", call. = FALSE)
    return(FALSE)
  }  
  lapply(wts, function(w){
    if (abs(length(w) - sum(w)) > 0.9){
      # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). 
      warning(sprintf("Empirical likelihood weights sum to %0.1f, which suggests the proposed null mean is on a face of the convex hull.", sum(w)), call. = FALSE)
    }
  })
  return(TRUE)
}


#' Eigenvalues divided to have squared sum 1 
#' Uses the fact that A*A squares the eigenvalues of A, and the trace of a matrix is the sum of the eigenvalues.
#' m A symmetric matric
#' @noRd
normL2evals <- function(m){
  ssq <- sum(diag(m %*% m))
  newm <- m / sqrt(ssq)
  newm <- makeSymmetric(newm) #remove computational inaccuracies
  return(newm)
}

#' @title Normalise so that Sum of Squared Eigenvalues is One
#' @description
#' Scales symmetric tensors so that the square of the eigenvalues sum to one.
#' @inheritParams test_multiplicity
#' @export
normalise_ss1 <- function(x){
  x <- as_fsm(x)
  if (ncol(x) == 6){#use fast method
    I2 <- x[, 1] * x[, 4] + x[, 4] * x[, 6] + x[,1] * x[,6] -
      x[, 2]^2 - x[, 5]^2 - x[,3]^2
    I1 <- x[, 1] + x[, 4] + x[, 6] #trace
    tr2 <- I1^2 - 2*I2
    return(as_fsm(x/sqrt(tr2)))
  }
  as_fsm(apply(x, 1, function(v){normL2evals(invvech(v))}, simplify = FALSE))
}

# ms is an TFORGE_fsm
normL2evals_sst <- normalise_ss1
