#' @title Test statistic for common eigenvalues across multiple samples
#' @description For the case that the eigenvectors are free, but the eigenvalues are fixed (single) sample.
#' @param mss List of samples
#' @details Corresponds to eqn 13 in the draft `tensors_4`.
stat_commonevals_ksample <- function(mss){
  esteval <- est_commonevals(mss)
  stat <- sum(vapply(mss, stat_specifiedevals, esteval, FUN.VALUE = 0.1))
  return(list(
    stat = stat,
    esteval = esteval
  ))
}

#' @title Estimate common eigenvalues from multiple samples
#' @details Eqn 24 in `tensors_4.pdf`, used for estimating eigenvalues for a multisample hypothesis test.
#' My implementation is surprisingly involved - there could be more elegant methods (or maybe my implementation is wrong).
est_commonevals <- function(mss){
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

