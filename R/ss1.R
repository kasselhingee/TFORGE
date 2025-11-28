#' @title Test for eigenvalues when sum of squared eigenvalues is 1
#' @description 
#' For a single sample of symmetric matrices where sum of squared eigenvalues  = 1, test eigenvalues of the population mean.
#' For multiple samples of symmetric matrices where sum of squared eigenvalues  = 1, test for equality of the eigenvalues of the population means.
#' The test statistic is calculated by `stat_ss1()`.
#' @details
#' Test hypotheses described below.
#' The sum of squared eigenvalues constraint forces the set of eigenvalues to lie on a sphere (or circle).
#' The test statistic accounts for this constraint by projecting eigenvalues onto a plane perpendicular to the direction of the sample average's eigenvalues.
#'
#' Weighted bootstrap calibration is used (see 'Weighted Bootstrapping' below).
#'
#' Eigenvalues must be distinct.
#' @inheritParams test_unconstrained
#' @inheritParams test_fixedtrace
#' @inheritSection test_unconstrained Hypotheses
#' @inheritSection test_fixedtrace Weighted Bootstrapping
#' @inherit test_unconstrained return
#' @references \insertAllCited{}
#' @export
test_ss1 <- function(x, evals = NULL, B = 1000, maxit = 25){
  x <- as_flat(x)
  stopifnot(has_ss1(x))
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  # check that evals satisfy ss1 constraint
  if (!is.null(evals)){if (abs(sqrt(sum(evals^2)) - 1) > sqrt(.Machine$double.eps)){stop("Square of evals do not sum to 1")}}
  
  if (has_fixedtrace(x)){warning("All tensors the same trace. Consider using test_ss1fixedtrace().")}
  
  if (B == "chisq"){
    df <- (dim_fsm_kfsm(x) - 1) * (length(x) - is.null(evals))
    return(chisq_calib(x, stat_ss1, df = df, evals = evals))
  }

  # statistic and null eigenvalues  
  t0 <- stat_ss1(x, evals = evals)
  d0 <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to d0)
  # also compute the bounds on free scalar c - See Eq37_cj_bound.pdf
  nullmeans <- lapply(x, elnullmean, d0 = d0, getcbound = TRUE)
  
  # maximise weights that empirical likelihood
  wts <- mapply(opt_el.test, ms = x, mu = nullmeans, maxit = maxit, SIMPLIFY = FALSE)

  
  #check the weights, if not good, exit
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
  
  res <- boot_calib(x, wts, 
                        stat = stat_ss1,
                        B = B,
                        evals = evals)
  return(res)
}

#' @rdname test_ss1
#' @export
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
  covars_unconstrained <- mapply(cov_evals_est, x = x, evecs = lapply(ess, "[[", "vectors"), av = mns, SIMPLIFY = FALSE)
  Deltas <- lapply(evalsav, function(d2){amaral2007Lemma1(d2/sqrt(sum(d2^2)))})
  Omega2s <- mapply(function(d2, covar, Delta){
    Delta %*% covar %*% t(Delta) / sum(d2^2)
  }, d2 = evalsav, covar = covars_unconstrained, Delta = Deltas, SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
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
    if (abs(sqrt(sum(evals^2)) - 1) > sqrt(.Machine$double.eps)){warning("Normalising eigenvalues supplied to stat_ss1() so that the squared eigenvalues sum to 1.")}
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
#' @description Compares whether the sum of the squared eigenvalues of the supplied symmetric matrices match each other using the property that the sum of the squared eigenvalues of `Y` equals the trace of `Y %*% Y`.
#' @inheritParams has_fixedtrace
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
# also compute the bounds on free scalar c - see Eq37_cj_bound.pdf
# @description The bounds on c are computed by exploiting the eigenvectors of the hypothesised mean:
# \deqn{cQdiag(\delta_0)Q^\top = \sum_{i=1}^n w_i Y_i}
# means that 
# \deqn{diag(c\delta_0) = \sum_{i=1}^n w_i Q^\top Y_i Q,}
# in particular
# \deqn{c\delta_0[1] = \sum_{i=1}^n w_i (Q^\top Y_i Q)[1,1]}
# and similar for 2nd, 3rd etc eigenvalue.
# Because the \eqn{w_i} create convex combinations, this means that
# \deqn{min((Q^\top Y_i Q) [1,1]) \leq c\delta_0[1] \leq max((Q^\top Y_i Q) [1,1]),}
# and similar for 2nd, 3rd etc eigenvalue.
# Currently this does not exploit the fact that \eqn{tr(Y_i Y_i) = 1} (and similar for \eqn{\delta_0}) and \eqn{c} is the Frobenius norm of \eqn{\sum_{i=1}^n w_i Y_i},
# which means that \eqn{c} must be between 0 and 1.
# @param ms a single sample of tensors
# @param d0 the NULL set of eigenvalues
# @param av The average of `ms`, if omitted then computed from `ms` directly (include to save computation time)
# @param evecs The eigenvectors of the average of `ms`. If omitted then computed from the average of `ms` directly (include to save computation time).
elnullmean <- function(ms, d0, av = NULL, evecs = NULL, getcbound = FALSE){
  if (is.null(av)){av <- mmean(ms)}
  if (is.null(evecs)){evecs <- eigen_desc(av)$vectors}
  nullmean <- evecs %*% diag(d0) %*% t(evecs)
  if (getcbound){
    # bounds for cj, which is the range of possible sqrt(sum of squared trace) of the resample means
    # A simpler method could be to just fix the bounds to [0,1].
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



#' For a single matrix, scale to have size 1 
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
#' @return A `TFORGE_fsm` object.
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

# small functino to help with writing unit tests
normalise_ss1_vec <- function(v){
  v / sqrt(sum(v^2))
}
