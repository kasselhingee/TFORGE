#' @noRd
#' @title Obsolete Tools for testing specified eigenvalues for a single sample
#' @description For the case that the eigenvectors are free and the eigenvalues are specified.
#' @param evals The eigenvalues are specified for single sample tests. In descending order.
#' @param ms Sample of matrices
NULL

#' @noRd
#' @title Test statistic (corresponds to eqn 13 in the draft `tensors_4`).
#' @param evecs Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of
stat_specifiedevals <- function(ms, evals, evecs = NULL){### OBSOLETE ###
  evals <- sort(evals, decreasing = TRUE)
  n <- nrow(ms)
  av <- mmean(as_fsm(ms))
  if (is.null(evecs)){
    av_eigenspace <- eigen_desc(av, symmetric = TRUE)
    d1 <- av_eigenspace$values
    evecs <- av_eigenspace$vectors
  } else {
    d1 <- diag(t(evecs) %*% av %*% evecs)
  }
  d0 <- evals
  V <- cov_evals_est(as_fsm(ms), evecs = evecs, av = av)
  out <- n * t(d1 - d0) %*% solve(V) %*% (d1 - d0)
  return(drop(out))
}


#' @noRd
#' @param B The number of bootstrap samples
test_specifiedevals <- function(ms, evals, B){
  ms <- as.mstorsst(ms)
  ms_std <- standardise_specifiedevals(ms, evals)
  res <- bootresampling(ms, ms_std, 
    stat = stat_specifiedevals,
    B = B,
    evals = evals)
  return(res)
}

