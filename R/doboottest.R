#' @title A function for helping perform single-sample bootstrap tests
#' @param ms List of matices - actual data
#' @param stdms List of matrices standardized to satisty the null.
#' @param stat Function to compute the statistic
#' @param B The number of bootstrap samples to use
singlesampletest <- function(ms, stdms, stat, B){
  t0 <- stat(ms)
  nullt <- replicate(B, stat(sample(stdms, replace = TRUE)))
  pval <- mean(nullt > t0)
  return(list(
    pval = pval,
    t0 = t0,
    nullt = nullt
  ))
}

