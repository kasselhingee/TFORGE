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
    evecs <- eigen_desc(av, symmetric = TRUE)$vectors
  }
  mcov <- mcovar(merr(ms, mean = av))

  # the V1 / V2 matrix for a single sample depending on whether evecs estimated or supplied
  indx <- expand.grid(1:nrow(evecs), 1:nrow(evecs))
  dupmat <- dup(nrow(evecs))
  vals <- mapply(cov_evals_inside, indx[,1], indx[,2],
                 MoreArgs = list(evecs = evecs, dupmat = dupmat, mcov = mcov))
  V <- matrix(NA, nrow = nrow(evecs), ncol = ncol(evecs))
  V[as.matrix(indx)] <- vals
  return(V)
}

cov_evals_inside <- function(j, k, evecs, dupmat, mcov){
  tmp <- evecs[, j] %*% t(evecs[, k])
  sum(diag(t(dupmat) %*% kronecker(tmp, tmp) %*% dupmat %*% mcov))
}



