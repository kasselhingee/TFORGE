#' @title Pivotal bootstrap test of mean eigenvalues
#' @description For a single sample of symmetric matrices, test eigenvalues of the population mean.
#' For multilpe samples of symmetric matrices, test for equality of the eigenvalues of the population means.
#' Eigenvalues must be distinct.
#' @details
#' For a single sample with `evecs` omitted, the null hypothesis is that the population mean has eigenvalues of `evals` without restriction on the eigenvectors. The alternative hypothesis here is that the eigenvalues are not `evals`.
#' If a `evecs` is supplied then both the null hypothesis and alternative hypothesis assume that the population mean has eigenvectors of `evecs`.
#'
#' For multiple samples, `evals` and `evecs` must be omitted. The null hypothesis is that the population means of each sample have the same eigenvalues.
#'
#' Bootstrap resampling is conducted from a population that satisfies the null hypothesis by transforming `x` with `standardise_specifiedevals()`.
#' The test statistic is calculated by `stat_unconstrained()`. 
#' @param x A single sample of symmetric matrices or multiple samples of symmetric matrices. See [`as_flat()`].
#' @param evals When `x` is a single sample, the null hypothesis is that the population mean has eigenvalues `evals`.
#' @param evecs For a single sample, specify eigenvectors to test under the assumption that the population mean's eigenvectors are the columns of `evecs`. The order of these eigenvectors matters and should be such that eigenvalues are in descending order.
#' @param B Number of bootstrap samples.
#' @return A `TFORGE` object with the eigenvalues of the null hypothesis in the `null_evals` attribute for `t0`. See [`bootresampling()`].
#' @examples
#' test_unconstrained(rsymm_norm(15, diag(c(3,2,1))), evals = c(3, 2, 1))
#' test_unconstrained(rsymm_norm(15, diag(c(3,2,1))), evals = c(3, 2, 1), evecs = diag(3))
#' test_unconstrained(list(rsymm_norm(15, diag(c(3,2,1))),
#'                         rsymm_norm(15, diag(c(3,2,1))))
#' @export
test_unconstrained <- function(x, evals = NULL, evecs = NULL, B = 1000){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since mss is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  if (!is.null(evecs) && (length(x) > 1)){stop("evecs specified for multisample not supported")}

  # user friendliness check for fixedtrace or ss1
  switch(4 - has_fixedtrace(x) - 2*has_ss1(x),
     warning("All tensors in x have the same trace and the sum of the squared eigenvalues is 1. Consider using test_ss1fixedtrace()"),
     warning("All tensors have a sum of squared eigenvalues of 1. Consider using test_ss1()"),
     warning("All tensors have the same trace. Consider using test_fixedtrace()"),
     NULL
     )
  
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

#' @rdname test_unconstrained
#' @export
stat_unconstrained <- function(x, evals = NULL, evecs = NULL){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  if (!is.null(evecs) && (length(x) > 1)){warning("evecs supplied to multisample, but stat_unconstrained() does not support providing evecs for a multisample")}

  avs <- lapply(x, mmean)
  ess <- lapply(avs, eigen_desc)

  #mean evals of samples
  if (is.null(evecs)){
    avevals <- lapply(ess, "[[", "values")
    avevecs <- lapply(ess, "[[", "vectors")
  } else {
    avevals <- lapply(avs, function(av){sort(diag(t(evecs) %*% av %*% evecs), decreasing = TRUE)})
    avevecs <- replicate(length(x), evecs, simplify = FALSE)
  }

  # covariance of evals (V matrices)
  Vs <- mapply(cov_evals_est, ms = x, evecs = avevecs, av = avs , SIMPLIFY = FALSE)
  
  # null evals
  if (is.null(evals)){
    d0 <- est_commonevals(x, Vs = Vs, evals = avevals)
    if (!all(order(d0, decreasing = TRUE) == 1:length(d0))){
      d0 <- descendingordererror(d0)
    }
  } else {
    d0 <- sort(evals, decreasing = TRUE)
  }

  # stat per sample
  persamplestat <- mapply(function(n, d1, V){
    n * t(d1 - d0) %*% solve_error(V) %*% (d1 - d0)
  },
  n = lapply(x, nrow),
  d1 = avevals,
  V = Vs,
  SIMPLIFY = FALSE)
  
  stat <- drop(purrr::reduce(persamplestat, `+`))
  attr(stat, "null_evals") <- drop(d0)
  return(stat)
}


#' Estimate eigenvalues in common to multiple sampler (Eqn 24 in `tensors_4.pdf`).
#' My implementation is surprisingly involved - there could be more elegant methods (or maybe my implementation is wrong).
#' Vs generated by applying cov_evals_est_est() to each sample
est_commonevals <- function(mss, Vs, evals){
  ns <- vapply(mss, nrow, FUN.VALUE = 1)
  Vs <- mapply(`/`, Vs, ns, SIMPLIFY = FALSE) #covariance of  Y_i eigenvalues to covariance of mean(Y) eigenvalues
  invVs <- lapply(Vs, solve_error)
  sum_invVs <- purrr::reduce(invVs, `+`)
  invVevals <- mapply(`%*%`, invVs, evals, SIMPLIFY = FALSE)
  sum_invVevals <- purrr::reduce(invVevals, `+`)
  return(drop(solve(sum_invVs) %*% sum_invVevals))
}

#' @rdname test_unconstrained
#' @export
standardise_specifiedevals <- function(ms, evals){
  ms <- as_flat(ms)
  evals <- sort(evals, decreasing = TRUE)
  av <- mmean(ms)
  errs <- merr(ms, mean = av)
  av_eigenspace <- eigen_desc(av, symmetric = TRUE)
  av_evecs <- av_eigenspace$vectors
  cen <- vech(av_evecs %*% diag(evals) %*% t(av_evecs))
  newms <- t(t(errs) + cen)
  class(newms) <- c("TFORGE_fsm", class(newms))
  return(newms)
}
