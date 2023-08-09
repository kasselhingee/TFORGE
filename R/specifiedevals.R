#' @name specifiedevals
#' @title Tools for testing specified eigenvalues for a single sample
#' @description For the case that the eigenvectors are free and the eigenvalues are specified.
#' @param evals The eigenvalues are specified for single sample tests. In descending order.
#' @param ms Sample of matrices
NULL

#' @describeIn specifiedevals Test statistic (corresponds to eqn 13 in the draft `tensors_4`).
#' @param evecs Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of
#' @export
stat_specifiedevals <- function(ms, evals, evecs = NULL){
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
  d0 <- evals
  V <- cov_evals(ms, evecs = evecs, av = av)
  out <- n * t(d1 - d0) %*% solve(V) %*% (d1 - d0)
  return(drop(out))
}



#' @describeIn specifiedevals Bootstrap test.
#' @param B The number of bootstrap samples
#' @export
test_specifiedevals <- function(ms, evals, B){
  ms_std <- standardise_specifiedevals(ms, evals)
  res <- singlesampletest(ms, ms_std, 
    stat = stat_specifiedevals,
    B = B,
    evals = evals)
  return(res)
}

#' @describeIn specifiedevals Standardise a sample to satisty the null hypothesis (i.e. the average has eigenvalues equal to `eval`).
#' @export
standardise_specifiedevals <- function(ms, evals){
  evals <- sort(evals, decreasing = TRUE)
  av <- mmean(ms)
  errs <- merr(ms, mean = av)
  av_eigenspace <- eigen(av, symmetric = TRUE)
  av_evecs <- av_eigenspace$vectors
  cen <- av_evecs %*% diag(evals) %*% t(av_evecs)
  newms <- lapply(errs, function(m) cen + m)
  return(newms)
}

#' General function for getting covariance of eigenvalues from a single sample
#' @param evecs If supplied these eigenvectors will be used instead of estimated eigenvectors. Each column is an eigenvector.
#' @param ms A single sample (list of matrices)
#' @param av To save computation, `av` can be passed if it is already computed.
#' @details Computes equation (11) of `tensors_4` with \eqn{C_0} replaced with the sample analogue. If `evecs` is not provided then the eigenvectors \eqn{q_{0i}} are replaced with the eigenvectors of the average of `ms`.
#' @return An estimated covariance matrix for the eigenvalues of `ms`.
#' @export
cov_evals <- function(ms, evecs = NULL, av = NULL){
  ms <- as.mstorsst(ms)
  stopifnot(inherits(ms, "sst"))
  if (is.null(av)){av <- mmean(ms)}
  if (is.null(evecs)){
    evecs <- eigen(av, symmetric = TRUE)$vectors
  }
  V <- cov_eval1_eval0(evecs, mcovar(merr(ms, mean = av)))
  return(V)
}

# the V1 / V2 matrix for a single sample depending on whether evecs estimated or supplied
cov_eval1_eval0 <- function(evecs, mcov){
  indx <- expand.grid(1:nrow(evecs), 1:nrow(evecs))
  dupmat <- dup(nrow(evecs))
  vals <- mapply(cov_eval1_eval0_inside, indx[,1], indx[,2],
    MoreArgs = list(evecs = evecs, dupmat = dupmat, mcov = mcov))
  out <- matrix(NA, nrow = nrow(evecs), ncol = ncol(evecs))
  out[as.matrix(indx)] <- vals
  return(out)
}

cov_eval1_eval0_inside <- function(j, k, evecs, dupmat, mcov){
  tmp <- evecs[, j] %*% t(evecs[, k])
  sum(diag(t(dupmat) %*% kronecker(tmp, tmp) %*% dupmat %*% mcov))
}


