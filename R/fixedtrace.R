# functions for testing eigenvalues when the trace of the matrices are fixed

#' @title Test whether the supplied sample(s) have fixed trace
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @details Credit to Hadley for much of this function: [https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector]
#' @export
#' @return `TRUE` or `FALSE`
hasfixedtrace <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as.mstorsst(x)
  if (inherits(x, "mst")){x <- unlist(x, recursive = FALSE)}
  traces <- vapply(x, function(y){sum(diag(y))}, FUN.VALUE = 1.64)
  tracerange <- range(traces) / mean(traces)
  isTRUE(all.equal(tracerange[1], tracerange[2], tolerance = tolerance))
}

#' @title Test statistic for given eigenvalues when trace is fixed.
#' @param evecs Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of
#' `t(evecs) %*% av %*% evecs`, where `av` is the average of `ms`.
#' @param ms A sample of matrices.
#' @param evals A set of hypothesised eigenvalues for the mean. `evals` must sum to the trace of the matrices.
#' @return A single numeric value.
stat_ss_fixedtrace <- function(ms, evals, evecs = NULL, NAonerror = TRUE){
  ms <- as.mstorsst(ms)
  stopifnot(hasfixedtrace((ms)))
  evals <- sort(evals, decreasing = TRUE)
  n <- length(ms)
  av <- mmean(ms)
  if (is.null(evecs)){
    av_eigenspace <- eigen(av, symmetric = TRUE)
    d1 <- av_eigenspace$values
    evecs <- av_eigenspace$vectors
  } else {
    d1 <- diag(t(evecs) %*% av %*% evecs)
  }
  if (!isTRUE(all.equal(sum(evals), sum(diag(av))))){stop("Provided evals do not sum to trace of average.")}
  d0 <- evals
  V <- cov_evals(ms, evecs = evecs, av = av)

  #project all eval covariances etc to plane orthogonal 1,1,1,1,1
  H <- helmertsub(ncol(av))
  Vproj <- H %*% V %*% t(H)
  
  erroraction <- function(e){
    if (!NAonerror){stop(e)}
    else {
      out <- NA * Vproj
      attr(out, "message") <- e$message
      return(out)
    }
  }
  Vprojinv <- tryCatch(solve(Vproj), error = erroraction)
  out <- n * t(d1 - d0) %*% t(H) %*% Vprojinv %*% H %*% (d1 - d0)
  return(drop(out))
}

#' @describeIn stat_ss_fixedtrace 
#' @param mss Multiple samples of matrices, all with the same trace. See [`as.mstorsst()`] for required structure.
stat_ms_fixedtrace <- function(mss, NAonerror = FALSE){
  mss <- as.mstorsst(mss)
  H <- helmertsub(ncol(mss[[1]][[1]]))
  ess <- lapply(mss, eigen)
  ns <- lapply(mss, length)
  
  #first get all eval precision matrices
  erroraction <- function(e){
    if (!NAonerror){stop(e)}
    else {
      out <- NA * diag(mss[[1]][[1]]$values)
      attr(out, "message") <- e$message
      return(out)
    }
  }
  precisions <- lapply(mss, function(ms){tryCatch(solve(H %*% cov_evals(ms) %*% t(H)), error = erroraction)})
  #get estimate of common evals
  sum_precisions <- purrr::reduce(precisions, `+`)
  precisionsbyevals <- mapply(`%*%`, precisions, lapply(ess, "[[", "values"), SIMPLIFY = FALSE)
  sum_precisionsbyevals <- purrr::reduce(precisionsbyevals, `+`)
  d0 <- drop(solve(sum_precisions) %*% sum_precisionsbyevals)
  
  #now compute the statistic, not using the single sample function to avoid redoing eigen()
  tmp <- mapply(function(n, d1, precision){
    n * t(d1 - d0) %*% t(H) %*% precision %*% H %*% (d1 - d0)
    },
    n = ns,
    d1 = lapply(ess, "[[", "values"),
    precision = precisions
  )
  stat <- purrr::reduce(tmp, `+`)
  attr(stat, "esteval") <- d0
  return(stat)
}
 
 
#' @describeIn stat_ss_fixedtrace Bootstrap test.
#' @param B The number of bootstrap samples
#' @export
test_ss_fixedtrace <- function(ms, evals, B, evecs = NULL){
  evals <- sort(evals, decreasing = TRUE)
  if (!isTRUE(all.equal(sum(evals), sum(diag(ms[[1]]))))){
    warning("Provided evals do not sum to trace of average.")
  }
  if (is.null(evecs)){
    av <- mmean(ms)
    es <- eigen(av, symmetric = TRUE)
    nullmean <- es$vectors %*% diag(evals) %*% t(es$vectors)
  } else {
    nullmean <- evecs %*% diag(evals) %*% t(evecs)
  }
  
  elres <- emplik::el.test(do.call(rbind, lapply(ms, vech)), vech(nullmean))
  if (abs(length(ms) - sum(elres$wts)) > 1E-2){
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is close or outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = Inf,
      nullt = NA,
      B = NA
    ))
  }
  
  res <- bootresampling(ms, elres$wts, 
                        stat = stat_ss_fixedtrace,
                        B = B,
                        evals = evals,
                        evecs = evecs)
  return(res)
}
  
