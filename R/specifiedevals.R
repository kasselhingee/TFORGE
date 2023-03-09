#' @title Test statistic for specified eigenvalues and single sample
#' @description For the case that the eigenvectors are free, but the eigenvalues are fixed (single) sample.
#' @param evals The eigenvalues are specified for single sample tests. In descending order.
#' @param ms Sample of matrices
#' @details Corresponds to eqn 13 in the draft `tensors_4`.
stat_commoneigenvals <- function(evals, ms){
  evals <- sort(evals, decreasing = TRUE)
  n <- length(ms)
  av <- mmean(ms)
  av_eigenspace <- eigen(av, symmetric = TRUE)
  av_evecs <- av_eigenspace$vectors
  d1 <- av_eigenspace$values
  d0 <- evals
  V <- cov_eval1_eval0(av_evecs, mcovar(merr(ms, mean = av)))
  out <- n * t(d1 - d0) %*% solve(V) %*% (d1 - d0)
  return(drop(out))
}

standardise_commoneigenvals <- function(evals, ms){
  evals <- sort(evals, decreasing = TRUE)
  av <- mmean(ms)
  errs <- merr(ms, mean = av)
  av_eigenspace <- eigen(av, symmetric = TRUE)
  av_evecs <- av_eigenspace$vectors
  cen <- av_evecs %*% diag(evals) %*% t(av_evecs)
  newms <- lapply(errs, function(m) cen + m)
  return(newms)
}

# the V1 / V2 matrix for a single sample
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


