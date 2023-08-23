#' @title Testing of unconstrained eigenvalues
#' @param x A multisample or single sample
#' @param evals For single sample, the NULL eigenvalues
#' @param NAonerror If TRUE, then failed matrix inversion result in NA, rather than an error.
#' @param evecs Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of `t(evecs) %*% av %*% evecs`, where `av` is the average of the sample.
#' @details The test statistic for a single sample corresponds to eqn 13 of draft `tensors_4.pdf`.
#' @return The statistic with the null eigenvalues used as an attribute `null_evals`.
#' @export
stat_unconstrained <- function(x, evals = NULL, evecs = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}

  avs <- lapply(x, mmean)
  ess <- lapply(avs, eigen)
  
  # null evals
  if (is.null(evals)){
    d0 <- est_commonevals(x, evecs = evecs, NAonerror = FALSE)
  } else {
    d0 <- sort(evals, decreasing = TRUE)
  }

  #mean evals of samples
  if (is.null(evecs)){
    avevals <- lapply(ess, "[[", "values")
    avevecs <- lapply(ess, "[[", "vectors")
  } else {
    avevals <- lapply(avs, function(av){sort(diag(t(evecs) %*% av %*% evecs), decreasing = TRUE)})
    avevecs <- replicate(length(x), evecs, simplify = FALSE)
  }


  # covariance of evals (V matrices)
  Vs <- mapply(cov_evals, ms = x, evecs = avevecs, av = avs , SIMPLIFY = FALSE)

  # stat per sample
  persamplestat <- mapply(function(n, d1, V){
    n * t(d1 - d0) %*% solve_NAonerror(V, NAonerror = NAonerror) %*% (d1 - d0)
  },
  n = lapply(x, length),
  d1 = avevals,
  V = Vs,
  SIMPLIFY = FALSE)
  
  stat <- drop(purrr::reduce(persamplestat, `+`))
  attr(stat, "null_evals") <- drop(d0)
  return(stat)
}


est_commonevals <- function(mss, evecs = NULL, NAonerror = FALSE){
  avs <- lapply(mss, mmean)
  if (is.null(evecs)){
    ess <- lapply(avs, eigen)
    evecs <- lapply(ess, "[[", "vectors")
    evals <- lapply(ess, "[[", "values")
  } else {
    evals <- lapply(avs, function(av){sort(diag(t(evecs) %*% av %*% evecs), decreasing = TRUE)})
    evecs <- replicate(length(mss), evecs, simplify = FALSE)
  }
  mcovars <- .mapply(function(a, b){
    mcovar(merr(a, mean = b))
  }, dots = list(a = mss, b = avs), MoreArgs = list())

  Vs <- .mapply(cov_eval1_eval0, 
          dots = list(evecs = evecs,
                      mcov = mcovars),
          MoreArgs = list())
  invVs <- lapply(Vs, solve_NAonerror, NAonerror = NAonerror)
  sum_invVs <- purrr::reduce(invVs, `+`)
  invVevals <- mapply(`%*%`, invVs, evals, SIMPLIFY = FALSE)
  sum_invVevals <- purrr::reduce(invVevals, `+`)
  return(drop(solve(sum_invVs) %*% sum_invVevals))
}
