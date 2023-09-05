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
#' @title Test for eigenvalues when trace is fixed.
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
#' @param evals If supplied the eigenvalues of the null hypothesis. When supplied `evals` must sum to the trace of the matrices. For the multisample statistic this should be `NULL` and the null evals estimated by the function.
stat_fixedtrace <- function(x, evals = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  stopifnot(hasfixedtrace(x))
  if (inherits(x, "sst")){mss <- as.mstorsst(list(x))}
  else {mss <- x}
  if (is.null(evals) && (length(mss) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(mss) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}

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
  precisions <- lapply(mss, function(ms){tryCatch(solve(H %*% cov_evals(ms) %*% t(H)), error = erroraction)})
  
  d1s <- lapply(ess, "[[", "values")
    rlang::warn("need to sort values in case of negative evals", .frequency = "once", .frequency_id = "devels")
  
  #get estimate of common evals for multisample situation
  if (is.null(evals)){
    sum_precisions <- purrr::reduce(precisions, `+`)
    precisionsbyevals <- mapply(function(A, B){A %*% H %*% B}, A = precisions, B = lapply(ess, "[[", "values"), SIMPLIFY = FALSE)
    rlang::warn("need to sort values in case of negative evals", .frequency = "once", .frequency_id = "devels")
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
 
 
#' @describeIn stat_fixedtrace Bootstrap test
#' @inheritParams stat_ss1fixedtrace
#' @export
test_fixedtrace <- function(x, evals = NULL, B, maxit = 25, sc = TRUE){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}

  t0 <- stat_fixedtrace(x, evals = evals, NAonerror = FALSE)
  estevals <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to estevals)
  nullmeans <- lapply(x, function(ms){
    av <- mmean(ms)
    evecs <- eigen(av)$vectors
    evecs %*% diag(estevals) %*% t(evecs)
  })
  
  # compute corresponding weights that lead to emp.lik.
  wts <- mapply(function(ms, nullmean){
    if (sc){
      scelres <- emplik(do.call(rbind, lapply(ms, vech)), vech(nullmean), itermax = maxit)
      if (!isTRUE(scelres$converged)){warning("emplik() did not converge")}
      wts <- as.vector(scelres$wts) * length(ms)
    } else {
      elres <- emplik::el.test(do.call(rbind, lapply(ms, vech)), vech(nullmean), maxit = maxit)
      if (elres$nits == maxit){warning(paste("Reached maximum iterations", maxit, "in el.test() at best null mean."))}
      wts <- elres$wts
    }
    wts
  }, ms = x, nullmean = nullmeans, SIMPLIFY = FALSE)

  #check the weights
  wtsums_discrepacies <- vapply(wts, function(x){abs(length(x) - sum(x))}, FUN.VALUE = 0.1)
  if (any(wtsums_discrepacies > 1E-2)){
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is close or outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    ))
  }
  
  res <- bootresampling(x, wts, 
                        stat = stat_fixedtrace,
                        B = B,
                        evals = evals)
  return(res)
}
