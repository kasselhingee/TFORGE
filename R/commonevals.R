#' @name commonevals
#' @title Tools for testing common eigenvalues across multiple samples
#' @description For the case that the eigenvalues are common across multiple samples (and unknown), and the eigenvectors are free.
#' __What is the appropriate way to handle singularities in resamples?__
#' @param mss List of samples. Each sample is itself a list of symmetric matrices.
NULL

#' @describeIn commonevals The test statistic (corresponds to eqn 13 in the draft `tensors_4`)
#' @param NAonerror If TRUE then when an error occurs NA values are returned. If error in estimating the eigenvalue then both return slots are `NA`. If an error occurs after eigenvalues are estimated, then just the `stat` slot is `NA`.
#' @export
stat_commonevals_ksample <- function(mss, NAonerror = FALSE){
  erroraction <- function(e){
    if (!NAonerror){stop(e)}
    else {
     out <- NA
     attr(out, "message") <- e$message
     return(out)
   }
  }
  esteval <- tryCatch(est_commonevals(mss), error = erroraction)
  stat <- tryCatch(sum(vapply(mss, stat_specifiedevals, esteval, FUN.VALUE = 0.1)), error = erroraction)
  return(list(
    stat = stat,
    esteval = esteval
  ))
}

#' @describeIn commonevals Estimate eigenvalues in common to multiple sampler (Eqn 24 in `tensors_4.pdf`).
#' My implementation is surprisingly involved - there could be more elegant methods (or maybe my implementation is wrong).
#' @export
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

#' @describeIn commonevals Bootstrap test of common eigenvalues
#' @param B Number of bootstrap samples
#' @export
test_commonevals <- function(mss, B){
  t0info <- stat_commonevals_ksample(mss)
  mss_std <- lapply(mss, standardise_specifiedevals, t0info$esteval)
  nullt <- replicate(B, { stat_commonevals_ksample(lapply(mss_std, sample, replace = TRUE), NAonerror = TRUE)$stat })
  if (any(is.na(nullt))){warning(sprintf("The statistic could not be calculated for %i bootstrap resamples.", sum(is.na(nullt))))}
  pval <- mean(nullt > t0info$stat, na.rm = TRUE)
  return(list(
   pval = pval,
   t0 = t0info$stat,
   nullt = nullt,
   esteval = t0info$esteval,
   B = B
  ))
}
