#' @title Test statistic for specified/common eigenvalues
#' @description For the case that the eigenvectors are free, but the eigenvalues are fixed (single) sample, or common.
#' @param evals The eigenvalues either specified for single sample tests, or estimated for multi-sample. In descending order.
#' @param ms Sample of matrices
#' @details Corresponds to eqn 13 in the draft `tensors_4`.

stat_commoneigenvals <- function(evals, ms){
  n <- length(ms)
  av <- mmean(ms)
  av_eigenspaces <- eigen(av, symmetric = TRUE)
  av_evecs <- av_eigenspaces$vectors
}

# the V1 / V2 matrix
cov_eval1_eval0 <- function(evecs, mcov){
  indx <- expand.grid(1:nrow(evecs), 1:nrow(evecs))
  dupmat <- dup(nrow(evecs))
  vals <- mapply(cov_eval1_eval0_inside, indx[,1], indx[,2],
    MoreArgs = list(evecs = evecs, dupmat = dupmat, mcov = mcov))
  out <- matrix(NA, nrow = nrow(evecs), ncol = ncol(evecs))
  out[indx] <- vals
  return(out)
}

cov_eval1_eval0_inside <- function(j, k, evecs, dupmat, mcov){
  tmp <- evecs[, j] %*% t(evecs[, k])
  sum(diag(t(dupmat) %*% kronecker(tmp, tmp) %*% dupmat %*% mcov))
}
