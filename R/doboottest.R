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

