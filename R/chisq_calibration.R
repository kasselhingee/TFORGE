#' Internal: Chi Squared Calibration for Testing
#' @description
#' Similar to [`bootresampling()`], but uses chi-squared calibration instead of bootstrapping.
#' @param df Degrees of freedom of the chi-squared distribution
#' @inheritParams bootresampling
#' @export
chisq_calib <- function(x, stat, df, ...){
  x <- as_flat(x)
  t0 <- stat(x, ...)
  pval <- 1-stats::pchisq(t0, df)
  attributes(pval) <- NULL
  
  out <- list(
    pval = pval,
    t0 = t0,
    df = df
  )
  class(out) <- c("TFORGE", class(out))
  return(out)
}
