#' Internal: Chi Squared Calibration for Testing
#' @description
#' Similar to `bootresampling()`, but uses a chi-squared distribution instead of an empirical distribution from bootstrapping.
#' @inheritParams bootresampling
#' @param df Degrees of freedom. If not supplied looks in `df` attribute in object returned by `stat`.
#' @export
chisq_calib <- function(x, stat, df = NULL, ...){
  x <- as_flat(x)
  t0 <- stat(x, ...)
  if (is.null(df)){df <- attr(t0, "df")}
  pval <- 1-stats::pchisq(t0, df = df)
  attributes(pval) <- NULL
  
  out <- list(
    pval = pval,
    t0 = t0,
    df = df
  )
  class(out) <- c("TFORGE", class(out))
  return(out)
}
