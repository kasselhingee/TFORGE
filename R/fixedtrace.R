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
#' @param evecs Optional argument for single sample statistic. Column vectors of eigenvalues, if supplied, the eigenvectors are considered fixed. In this case the \eqn{\delta_1} in the statistic is the diagonal of
#' `t(evecs) %*% av %*% evecs`, where `av` is the average of `ms`.



#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
#' @param evals If supplied the eigenvalues of the null hypothesis and `evals` must sum to the trace of the matrices. For the multisample statistic this should be `NULL` and is estimated within the function.
stat_fixedtrace <- function(x, evals = NULL, evecs = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  stopifnot(hasfixedtrace(x))
  if (is.null(evals) && inherits(x, "sst")){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && inherits(x, "mst")){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  if (!is.null(evecs) && inherits(x, "mst")){warning("evecs supplied for multisample situation supplied. This is unusual.")}
  if (inherits(x, "sst")){mss <- as.mstorsst(list(x))}
  else {mss <- x}
  H <- helmertsub(ncol(mss[[1]][[1]]))
  ess <- lapply(mss, function(ms){eigen(mmean(ms))})
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
  precisions <- lapply(mss, function(ms){tryCatch(solve(H %*% cov_evals(ms, evecs = evecs) %*% t(H)), error = erroraction)})
  
  # if evecs supplied the as per tensors_4.pdf T_1 statistic, using supplied evecs to test eigenvalues
  if (is.null(evecs)){
    d1s <- lapply(ess, "[[", "values")
  } else {
    d1s <- lapply(mss, function(ms){diag(t(evecs) %*% mmean(ms) %*% evecs)})
  }
  
  #get estimate of common evals for multisample situation
  if (is.null(evals)){
    sum_precisions <- purrr::reduce(precisions, `+`)
    precisionsbyevals <- mapply(function(A, B){A %*% H %*% B}, A = precisions, B = lapply(ess, "[[", "values"), SIMPLIFY = FALSE)
    sum_precisionsbyevals <- purrr::reduce(precisionsbyevals, `+`)
    d0proj <- drop(solve(sum_precisions) %*% sum_precisionsbyevals)
    d0 <- (t(H) %*% d0proj) + mean(diag(mss[[1]][[1]])) #convert projected evals back to p-dimensions, then shift to give correct trace.
    if (!all(order(d0) == length(d0):1)){
      d0 <- sort(d0)
      warning("Estimated common eigenvalues are not in descending order and has been reordered.")
    }
  } else {
    if (!isTRUE(all.equal(sum(evals), sum(diag(mss[[1]][[1]]))))){stop("Provided evals do not sum to trace of observations.")}
    d0 <- sort(evals, decreasing = TRUE)
  }
  
  #now compute the statistic.
  tmp <- mapply(function(n, d1, precision){
    n * t(d1 - d0) %*% t(H) %*% precision %*% H %*% (d1 - d0)
    },
    n = ns,
    d1 = d1s,
    precision = precisions,
    SIMPLIFY = FALSE
  )
  stat <- drop(purrr::reduce(tmp, `+`))
  
  attr(stat, "null_evals") <- drop(d0)
  if (is.null(evals)){attr(stat, "null_evals_proj") <- drop(d0proj)}
  return(stat)
}
 
 
#' @describeIn stat_fixedtrace Bootstrap test.
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
                        stat = stat_fixedtrace,
                        B = B,
                        evals = evals,
                        evecs = evecs)
  return(res)
}

#' @describeIn stat_fixedtrace Bootstrap test.
#' @export
test_ms_fixedtrace <- function(mss, B){
  mss <- as.mstorsst(mss)
  t0 <- stat_fixedtrace(mss, NAonerror = FALSE)
  evals <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to evals)
  nullmeans <- lapply(mss, function(ms){
    av <- mmean(ms)
    evecs <- eigen(av)$vectors
    evecs %*% diag(evals) %*% t(evecs)
  })
  
  # compute corresponding weights that lead to emp.lik.
  wts <- mapply(function(ms, nullmean){
    elres <- emplik::el.test(do.call(rbind, lapply(ms, vech)), vech(nullmean))
    elres$wts
  }, ms = mss, nullmean = nullmeans, SIMPLIFY = FALSE)

  #check the weights
  wtsums_discrepacies <- vapply(wts, function(x){abs(length(x) - sum(x))}, FUN.VALUE = 0.1)
  if (any(wtsums_discrepacies > 1E-2)){
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is close or outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = Inf,
      nullt = NA,
      B = NA
    ))
  }
  
  res <- bootresampling(mss, wts, 
                        stat = stat_fixedtrace,
                        B = B)
  return(res)
}
