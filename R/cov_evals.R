#' General function for computing covariance of eigenvalues.
#' @param evecs Matrix of eigenvectors of the mean as columns.
#' @param mcov Covariance of `vech(Y)`, where `Y` is the random matrix.
#' @details Computes equation (11) of `tensors_4` with \eqn{C_0}.
#' @useDynLib TFORGE, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @export
cov_evals <- function(evecs, mcov){
  # the V1 / V2 matrix for a single sample depending on whether evecs estimated or supplied
  indx <- rbind(t(combn(1:nrow(evecs), 2)), #this avoids repeating elements that are symmetric
                cbind(1:nrow(evecs), 1:nrow(evecs)))
  dupmat <- dup(nrow(evecs)) # is sparse so could be even faster
  vals <- mapply(cov_evals_inside_cpp, 
                 lapply(indx[,1], function(i) {evecs[, i]}),
                 lapply(indx[,2], function(i) {evecs[, i]}),
                 MoreArgs = list(dupmat = dupmat, mcov = mcov))
  V <- matrix(NA, nrow = nrow(evecs), ncol = ncol(evecs))
  V[as.matrix(indx)] <- vals
  V[lower.tri(V)] <- t(V)[lower.tri(V)]
  return(V)
}


# @param ms A single sample (list of matrices)
# @param av To save computation, `av` can be passed if it is already computed.
# @details If `evecs` is not provided then the eigenvectors \eqn{q_{0i}} are replaced with the eigenvectors of the average of `ms`.
# @return An estimated covariance matrix for the eigenvalues of `ms`.
cov_evals_est <- function(ms, evecs = NULL, av = NULL){
  ms <- as.mstorsst(ms)
  stopifnot(inherits(ms, "sst"))
  if (is.null(av)){av <- mmean(ms)}
  if (is.null(evecs)){
    evecs <- eigen_desc(av, symmetric = TRUE)$vectors
  }
  mcov <- mcovar(merr(ms, mean = av))
  if (all(mcov == 0)){
    stop(structure(
      class = c("zerocovariance", "error", "condition"),
      list(message = "covariance is zero",
           call =  sys.call(0))
    ))
  }

  cov_evals(evecs, mcov)
}

cov_evals_inside <- function(vecj, veck, dupmat, mcov){
  tmp <- vecj %*% t(veck)
  sum(diag(t(dupmat) %*% kronecker(tmp, tmp) %*% dupmat %*% mcov))
}


