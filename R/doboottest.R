#' @title A function for helping perform single-sample and k-sample bootstrap tests
#' @param x Observations as a list of matices, or list of list of matrices.
#' @param stdx List of matrices standardized to satisfy the null OR sampling weights for each matrix in `x`, in the same structure as `x`.
#' @param stat Function to compute the statistic
#' @param B The number of bootstrap samples to use
#' @param ... Passed to `stat`
#' @export
bootresampling <- function(x, stdx, stat, B, ...){
  x <- as.mstorsst(x)
  t0 <- stat(x, ...)
  exargs <- list(...)
  if (inherits(x, "mst")){
    if (inherits(stdx[[1]][[1]], "numeric")){
      #stdx is weights for an mst because first element of first sample is not a matrix/array, but just a numeric
      nullt <- replicate(B, do.call(stat, c(list(multisample(x, prob = stdx)), exargs)))
    } else {
      nullt <- replicate(B, do.call(stat, c(list(multisample(stdx)), exargs)))
    }
  } else if (inherits(x, "sst")){
    if (inherits(stdx[[1]], "numeric")){
      #stdx is weights for an sst because first element sample is not a matrix/array, but just a numeric
      nullt <- replicate(B, do.call(stat, c(list(sample(x, prob = stdx, replace = TRUE)), exargs)))
    } else {
      nullt <- replicate(B, do.call(stat, c(list(sample(stdx, replace = TRUE)), exargs)))
    }
  }
  
  if (any(is.na(nullt))){warning(sprintf("The statistic could not be calculated for %i bootstrap resamples.", sum(is.na(nullt))))}
  
  pval <- mean(nullt > t0)
  return(list(
    pval = pval,
    t0 = t0,
    nullt = nullt,
    B = B
  ))
}

#equivalent of sample, but for multiple samples
#' @param x a `mst`
#' @param w weights. If present, must have the same structure as `x`
multisample <- function(x, prob = NULL){
  if (is.null(prob)){return(lapply(x, sample, replace = TRUE))}
  else {
    stopifnot(length(prob) == length(x))
    stopifnot(all(vapply(x, length, 2) == vapply(prob, length, 2)))
    return(mapply(sample, x, prob = prob, MoreArgs = list(replace = TRUE)))
  }
}


