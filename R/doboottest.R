#' @title A function for helping perform single-sample bootstrap tests
#' @param ms List of matices - actual data
#' @param stdms List of matrices standardized to satisty the null.
#' @param stat Function to compute the statistic
#' @param B The number of bootstrap samples to use
#' @param ... Passed to `stat`
singlesampletest <- function(ms, stdms, stat, B, ...){
  t0 <- stat(ms, ...)
  exargs <- list(...)
  nullt <- replicate(B, do.call(stat, c(list(ms = sample(stdms, replace = TRUE)),
                                        exargs)))
  pval <- mean(nullt > t0)
  return(list(
    pval = pval,
    t0 = t0,
    nullt = nullt,
    B = B
  ))
}

#' @title A function for helping perform single-sample and k-sample bootstrap tests
#' @param x Observations as a list of matices, or list of list of matrices.
#' @param stdx List of matrices standardized to satisty the null OR sampling weights for each matrix in `x`, in the same structure as `x`.
#' @param stat Function to compute the statistic
#' @param B The number of bootstrap samples to use
#' @param ... Passed to `stat`
bootresampling <- function(x, stdx, stat, B, ...){
  x <- as.mstorsst(x)
  t0 <- stat(x, ...)
  exargs <- list(...)
  nullt <- replicate(B, do.call(stat, c(list(ms = sample(stdms, replace = TRUE)),
                                        exargs)))
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
multisample <- function(x, w = NULL){
  if (is.null(w)){return(lapply(x, sample, replace = TRUE))}
  else {
    stopifnot(length(w) == length(x))
    stopifnot(all(vapply(x, length, 2) == vapply(w, length, 2)))
    return(mapply(sample, x, prob = w, MoreArgs = list(replace = TRUE)))
  }
}


