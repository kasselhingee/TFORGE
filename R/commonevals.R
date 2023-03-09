#' @title Test statistic for specified/common eigenvalues
#' @description For the case that the eigenvectors are free, but the eigenvalues are fixed (single) sample, or common.
#' @param evals The eigenvalues either specified for single sample tests, or estimated for multi-sample. In descending order.
#' @param ms Sample of matrices
#' @details Corresponds to eqn 13 in the draft `tensors_4`.

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

#' @title The test statistic for a set of samples
stat_commonevals_ksample <- function(mss){
  esteval <- est_commoneigenvals(mss)
  stat <- sum(vapply(mss, function(Y){stat_commoneigenvals(esteval, Y)}, FUN.VALUE = 0.1))
  return(list(
    stat = stat,
    esteval = esteval
  ))
}

#' @title Estimate common eigenvalues from multiple samples
#' @details Eqn 24 in `tensors_4.pdf`, used for estimating eigenvalues for a multisample hypothesis test.
#' My implementation is surprisingly involved - there could be more elegant methods (or maybe my implementation is wrong).
#' @param mss A list of samples
est_commoneigenvals <- function(mss){
  avs <- lapply(mss, mmean)
  ess <- lapply(avs, eigen)
  mcovars <- .mapply(function(a, b){
    mcovar(merr(a, mean = b))
  }, dots = list(a = mss, b = avs), MoreArgs = list())

  Vs <- .mapply(cov_eval1_eval0, 
          dots = list(evecs = lapply(ess, "[[", "vectors"),
                      mcov = mcovars),
          MoreArgs = list())
  invVs <- lapply(Vs, solve)
  sum_invVs <- purrr::reduce(invVs, `+`)
  invVevals <- mapply(`%*%`, invVs, lapply(ess, "[[", "values"), SIMPLIFY = FALSE)
  sum_invVevals <- purrr::reduce(invVevals, `+`)
  return(drop(solve(sum_invVs) %*% sum_invVevals))
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


