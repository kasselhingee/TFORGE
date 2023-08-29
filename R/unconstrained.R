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
  if (!is.null(evecs) && (length(x) > 1)){warning("evecs supplied to multisample, but stat_unconstrained() does not support providing evecs for a multisample")}

  avs <- lapply(x, mmean)
  ess <- lapply(avs, eigen)
  
  # null evals
  if (is.null(evals)){
    d0 <- est_commonevals(x, evecs = evecs, NAonerror = NAonerror)
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

#' @describeIn stat_unconstrained Bootstrap test using `stat_unconstrained()`.
#' @param B Number of bootstrap samples.
#' @export
test_unconstrained <- function(x, evals = NULL, evecs = NULL, B){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since mss is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  if (!is.null(evecs) && (length(x) > 1)){stop("evecs specified for multisample not supported")}
  
  if (is.null(evals)){#estimate common evals using stat_unconstrained()
    t0info <- stat_unconstrained(x, evecs = evecs)
    estevals <- attr(t0info, "null_evals") #estevals name here because don't pass estimated evals to bootresampling
  } else {
    estevals <- evals
  }
  x_std <- lapply(x, standardise_specifiedevals, estevals)
  
  res <- bootresampling(x, x_std, 
                        stat = stat_unconstrained,
                        B = B,
                        evals = evals,
                        evecs = evecs)
  return(res)
}

#' @describeIn stat_unconstrained Estimate eigenvalues in common to multiple sampler (Eqn 24 in `tensors_4.pdf`).
#' My implementation is surprisingly involved - there could be more elegant methods (or maybe my implementation is wrong).
#' @export
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

#' @describeIn stat_unconstrained Standardise a sample to satisfy the null hypothesis (i.e. the average has eigenvalues equal to `eval`).
#' @export
standardise_specifiedevals <- function(ms, evals){
  ms <- as.mstorsst(ms)
  evals <- sort(evals, decreasing = TRUE)
  av <- mmean(ms)
  errs <- merr(ms, mean = av)
  av_eigenspace <- eigen(av, symmetric = TRUE)
  av_evecs <- av_eigenspace$vectors
  cen <- av_evecs %*% diag(evals) %*% t(av_evecs)
  newms <- lapply(errs, function(m) cen + m)
  newms <- lapply(newms, function(out){out[lower.tri(out)] <- out[upper.tri(out)]; out}) #to remove machine differences
  class(newms) <- c(class(newms), "sst")
  return(newms)
}
