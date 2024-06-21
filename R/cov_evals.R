#' Eigenvalue covariance from eigenvectors and element-wise covariance
#' @description
#' For a random symmetric matrix `Y`, the covariance of the eigenvalues of `Y` is calculated from the covariance of the elements and the eigenvectors of the mean of `Y`.
#' This calculation is valid when the eigenvalues are distinct.
#' @param evecs Matrix of eigenvectors of the mean of `Y`. Eigenvectors as columns and ordered according to descending eigenvalue.
#' @param mcov Covariance of `vech(Y)`, where `Y` is the random matrix.
#' @details Uses `RcppArmadillo` for a slight increase in speed.
# Computes equation (11) of `tensors_4` with \eqn{C_0}.
#' @useDynLib TFORGE, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @return A covariance matrix corresponding to the eigenvalues of the mean of `Y` in descending order.
#' @export
cov_evals <- function(evecs, mcov){
  # the V1 / V2 matrix for a single sample depending on whether evecs estimated or supplied
  indx <- rbind(t(utils::combn(1:nrow(evecs), 2)), #this avoids repeating elements that are symmetric
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


# @param x A single sample (list of matrices)
# @param av To save computation, `av` can be passed if it is already computed.
# @details If `evecs` is not provided then the eigenvectors \eqn{q_{0i}} are replaced with the eigenvectors of the average of `x`.
# @return An estimated covariance matrix for the eigenvalues of `x`.
cov_evals_est <- function(x, evecs = NULL, av = NULL){
  x <- as_flat(x)
  stopifnot(inherits(x, "TFORGE_fsm"))
  if (is.null(av)){av <- mmean(x)}
  if (is.null(evecs)){
    evecs <- eigen_desc(av, symmetric = TRUE)$vectors
  }
  mcov <- mcovar(merr(x, mean = av))
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


