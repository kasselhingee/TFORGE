# functions for testing eigenvalues when the trace of the matrices are fixed

#' @title Test whether the supplied sample(s) have fixed trace
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @details Credit to Hadley for much of this function: [https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector]
#' @export
#' @return `TRUE` or `FALSE`
hasfixedtrace <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as.mstorsst(x)
  if (inherits(x, "mst")){x <- unlist(x, recursive = FALSE)}
  traces <- vapply(x, function(y){sum(diag(y))}, FUN.VALUE = 1.64)
  tracerange <- range(traces) / mean(traces)
  isTRUE(all.equal(tracerange[1], tracerange[2], tolerance = tolerance))
}

#' @title Test statistic for given eigenvalues when trace is fixed.
#' @param evecs Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of
#' `t(evecs) %*% av %*% evecs`, where `av` is the average of `ms`.
#' @param ms A sample of matrices.
#' @param evals A set of hypothesised eigenvalues for the mean. `evals` must sum to the trace of the matrices.
#' @return A single numeric value.
stat_ss_fixedtrace <- function(ms, evals, evecs = NULL){
  browser()
  ms <- as.mstorsst(ms)
  stopifnot(hasfixedtrace((ms)))
  evals <- sort(evals, decreasing = TRUE)
  n <- length(ms)
  av <- mmean(ms)
  if (is.null(evecs)){
    av_eigenspace <- eigen(av, symmetric = TRUE)
    d1 <- av_eigenspace$values
    evecs <- av_eigenspace$vectors
  } else {
    d1 <- diag(t(evecs) %*% av %*% evecs)
  }
  if (!isTRUE(all.equal(sum(evals), sum(diag(av))))){stop("Provided evals do not sum to trace of average.")}
  d0 <- evals
  V <- cov_evals(ms, evecs = evecs, av = av)

  #project all eval covariances etc to plane orthogonal 1,1,1,1,1
  H <- helmertsub(ncol(av))
  Vproj <- H %*% V %*% t(H)
  out <- n * t(d1 - d0) %*% t(H) %*% solve(Vproj) %*% H %*% (d1 - d0)
  return(drop(out))
}
  
  
