#' @title Obsolete Tools for testing common eigenvalues across multiple samples
#' @description For the case that the eigenvalues are common across multiple samples (and unknown), and the eigenvectors are free.
#' __What is the appropriate way to handle singularities in resamples?__
#' @param mss List of samples. Each sample is itself a list of symmetric matrices.

#' @describeIn commonevals The test statistic (corresponds to eqn 13 in the draft `tensors_4`)
#' @param NAonerror If TRUE then when an error occurs NA values are returned. If error is before estimating the eigenvalues then the attribute `est_eval` is also `NA`.
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
  attr(stat, "esteval") <- esteval
  return(stat)
}

test_commonevals <- function(mss, B){
  t0info <- stat_commonevals_ksample(mss)
  mss_std <- lapply(mss, standardise_specifiedevals, attr(t0info, "esteval"))
  out <- bootresampling(mss, mss_std, stat_commonevals_ksample, B = B, NAonerror = TRUE)
  return(c(out, list(esteval = attr(t0info, "esteval"))))
}


