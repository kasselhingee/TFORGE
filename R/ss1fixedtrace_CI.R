#' @title Confidence interval for 3x3 Tensors with trace 0 and sum of sq eigenvalues equal to 1
#' @param x A sample of 3x3 tensors.
#' @param alpha Significance level
#' @param B Number of bootstrap resamples.
#' @return A list:
#' + `est`: point estimate of eigenvalues of the mean tensor
#' + `lower`: lower end of the confidence interval
#' + `upper`: upper end of the confidence interval
#' + `Omega`: The estimated covariance of the (projected) eigenvalues
#' + `threshold`: The threshold on the statistic used.
#' @export
conf_ss1fixedtrace <- function(x, alpha = 0.05, B = 1000){
  # checks
  x <- as.sst(x)
  stopifnot(ncol(x) == 6)
  stopifnot(hasfixedtrace(x))
  stopifnot(hasss1(x))
  
  # get sample mean
  av <- mmean(x)
  es <- eigen_desc(av)
  naveval <- es$values/sqrt(sum(es$values^2))
  
  # resample and get statistic each time
  res <- bootresampling(x, x, stat = stat_ss1fixedtrace, B = B, evals = naveval)
  warning("check: need to ignore original sample")
  # get 1-alpha quantile of the resampled statistics
  statthreshold <- quantile(res$nullt, probs = 1-alpha, names = FALSE, type = 1)
  
  # get V from x
  covar_unconstrained <- cov_evals_est(ms = x, evecs = es$vectors, av = av)
  
  
  # get lower and upper bound by solving quadratic formula
  # a building block
  A0 <- matrix(c( 0, 1,-1,
                  -1, 0, 1,
                  1,-1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  # the x^2 in my notes
  exsq <- statthreshold * t(naveval) %*% t(A0) %*% covar_unconstrained %*% A0 %*% naveval/9
  ex_upper <- sqrt(drop(exsq))
  
  # check bounds do not include multiplicity situations - because ex interval is symmetric about 0, sign not needed
  eval_21mult <- c(1, 1, -2)/sqrt(6)
  ex_21mult <- drop(tan(acos(t(eval_21mult) %*% naveval))/sqrt(3)) #sqrt(3) here is the size of A0 %*% naveval
  #checked: drop(-ex_21mult) * drop(A0 %*% naveval) + naveval has 2-1 multiplicity
  intersect_21mult <- abs(ex_21mult) <= ex_upper
  
  eval_12mult <- c(2, -1, -1)/sqrt(6)
  ex_12mult <- drop(tan(acos(t(eval_12mult) %*% naveval))/sqrt(3)) #sqrt(3) here is the size of A0 %*% naveval
  #checked: ex_12mult * drop(A0 %*% naveval) + naveval #has 1 - 2 multiplicity
  intersect_12mult <- abs(ex_12mult) <= ex_upper

  if (intersect_21mult | intersect_12mult){
    desc <- switch(intersect_21mult + 2 * intersect_12mult,
                   "largest two eigenvalues",
                   "smallest two eigenvalues",
                   "largest two eigenvalues or smallest two eigenvalues")
    warning(sprintf("Confidence region includes eigenvalues where the %s are not in descending order.", desc))
  } 
  
  # warn if interval is getting close to half a circle (where the projection A0 %*% naveval shouldn't work so well)
  if (drop(tan(pi/3)/sqrt(3)) <= ex_upper){
    warning("Confidence interval is very long (more than 1/3 of a circle) - projection onto the tangent at eigenvalue of mean may be poor.")
  }
  
  # return bounds as full vectors
  lower <- -ex_upper * A0 %*% naveval + naveval
  lower <- drop(lower / sqrt(sum(lower^2)))
  upper <- ex_upper * A0 %*% naveval + naveval
  upper <- drop(upper / sqrt(sum(upper^2)))
  
  return(list(lower = lower, upper = upper))
}
  

# and a test function - see if lower and upper points project to either side of the tangent at the test evals
#' @describeIn conf_ss1fixedtrace Test whether a particular set of eigenvalues `evals` lies in the confidence region returned by [conf_ss1fixedtrace()].
#' @param evals A set of eigenvalues with trace of zero and sum of squares of 1.
#' @param cr A confidence region returned by [conf_ss1fixedtrace()].
#' @export
conf_ss1fixedtrace_inregion <- function(evals, cr){
  stopifnot(sum(evals) <= sqrt(.Machine$double.eps))
  evals <- evals/sqrt(sum(evals^2))
  A0 <- matrix(c( 0, 1,-1,
                  -1, 0, 1,
                  1,-1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  lower_proj <- drop(t(A0 %*% evals) %*% cr$lower)
  upper_proj <- drop(t(A0 %*% evals) %*% cr$upper)
  sameside <- (drop(evals %*% cr$lower) > 0) & (drop(evals %*% cr$upper) > 0)
  return(sameside & (lower_proj <= sqrt(.Machine$double.eps)) & (upper_proj >= -sqrt(.Machine$double.eps)))
}
