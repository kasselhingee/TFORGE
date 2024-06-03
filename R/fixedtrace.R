# functions for testing eigenvalues when the trace of the matrices are fixed

#' @title Test whether the supplied sample(s) have fixed trace
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @details Credit to Hadley for much of this function: [https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector]
#' @export
#' @return `TRUE` or `FALSE`
hasfixedtrace <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as.mstorsst(x)
  if (inherits(x, "mst")){x <- do.call(rbind, x)}
  diagels <- isondiag_vech(x[1, ])
  traces <- rowSums(x[, diagels])
  tracerange <- range(traces)
  isTRUE(all.equal(tracerange[1], tracerange[2], tolerance = tolerance))
}

#' @title Test for eigenvalues when trace is fixed.
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
#' @param evals If supplied the eigenvalues of the null hypothesis. When supplied `evals` must sum to the trace of the matrices. For the multisample statistic this should be `NULL` and the null evals estimated by the function.
#' @export
stat_fixedtrace <- function(x, evals = NULL){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){mss <- as.mstorsst(list(x))}
  else {mss <- x}
  if (is.null(evals) && (length(mss) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(mss) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}

  H <- helmertsub(dimfromvech(mss[[1]][1,]))
  avs <- lapply(mss, mmean)
  ess <- lapply(avs, function(av){eigen_desc(av)})
  ns <- lapply(mss, nrow)
  
  #first get all eval precision matrices 
  ## below could be faster by passing evecs to cov_evals_est()
  ## stop using solve_error, or atleast the parameter
  precisions <- mapply(function(ms, evecs, av){solve_error(cov_evals_ft(ms, H = H, evecs = evecs, av = av))},
                       ms = mss,
                       evecs = lapply(ess, "[[", "vectors"),
                       av = avs,
                       SIMPLIFY = FALSE)
  d1s <- lapply(ess, "[[", "values")
  
  #get estimate of common evals for multisample situation
  if (is.null(evals)){
    precisions_mean <- mapply(`*`, ns, precisions, SIMPLIFY = FALSE)#precisions of the means, not the Y
    sum_precisions <- purrr::reduce(precisions_mean, `+`)
    precisionsbyevals <- mapply(function(A, B){A %*% H %*% B}, A = precisions_mean, B = lapply(ess, "[[", "values"), SIMPLIFY = FALSE)
    sum_precisionsbyevals <- purrr::reduce(precisionsbyevals, `+`)
    d0proj <- drop(solve_error(sum_precisions) %*% sum_precisionsbyevals)
    d0 <- (t(H) %*% d0proj) + mean(diag(invvech(mss[[1]][1, ]))) #convert projected evals back to p-dimensions, then shift to give correct trace.
    if (!all(order(d0, decreasing = TRUE) == 1:length(d0))){
      d0 <- descendingordererror(d0)
    }
  } else {
    if (!isTRUE(all.equal(sum(evals), sum(diag(invvech(mss[[1]][1, ])))))){stop("Provided evals do not sum to trace of observations.")}
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
test_fixedtrace <- function(x, evals = NULL, B, maxit = 25){
  x <- as.mstorsst(x)
  stopifnot(hasfixedtrace(x))
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  if (!is.null(evals)){
    if (all(abs(evals - evals[1]) < sqrt(.Machine$double.eps))){
      warning("Supplied evals are equal to each other so test is testing for isotropy. test_multiplicity() usually has better behaviour for testing isotropy.")
    }
  }


  t0 <- stat_fixedtrace(x, evals = evals)
  estevals <- attr(t0, "null_evals")
  
  # compute means that satisfy the NULL hypothesis (eigenvalues equal to estevals)
  nullmeans <- lapply(x, function(ms){
    av <- mmean(ms)
    evecs <- eigen_desc(av)$vectors
    evecs %*% diag(estevals) %*% t(evecs)
  })
  
  # compute corresponding weights that lead to emp.lik.
  wts <- mapply(function(ms, nullmean){
    scelres <- emplik(ms, vech(nullmean), itermax = maxit)
    if (!isTRUE(scelres$converged)){warning("emplik() did not converge, which usually means that the proposed null mean is outside the convex hull of the data")}
    wts <- as.vector(scelres$wts) * nrow(ms)
    wts
  }, ms = x, nullmean = nullmeans, SIMPLIFY = FALSE)

  #check the weights
  if (!wtsokay(wts)){
    out <- list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    )
    class(out) <- c("tensorboot", class(out))
    return(out)
  }

  res <- bootresampling(x, wts, 
                        stat = stat_fixedtrace,
                        B = B,
                        evals = evals)
  return(res)
}

#' Diagonal elements projected onto hyperplane through origin orthogonal to (1,1,..., 1)
#' m A symmetric matrix
projtrace <- function(m){ #project to have trace 0
  diags <- diag(m)
  ones <- rep(1, length(diags))/sqrt(length(diags))
  newdiag <- diags - drop(diags %*% ones) * ones
  diag(m) <- newdiag
  return(m)
}
projtrace_sst <- function(ms){
  diagels <- isondiag_vech(ms[1, ])
  H <- helmert(sum(diagels))
  projmat <- t(H) %*% diag(c(0,rep(1, sum(diagels)-1))) %*% H
  diags <- ms[, diagels, drop = FALSE]
  ms[, diagels] <- diags %*% t(projmat)
  return(ms)
}

#' @title Scale tensors to have trace of one.
#' @description Scales tensors by their trace, so that resulting tensors have a trace of one.
#' The method will create `Inf` values for tensors that have a trace of zero.
# This is how magnetic susceptibility tensors are scaled to examine isotropy in paleomagnetism (Tauxe, 2010 Essentials of Paleomagnetism, p248).
# I think will usually break the symmetry of the eigenvalues about their mean.
#' @param x A symmetric tensor, or a set of symmetric tensors as an `sst`.
#' @export
normalise_trace <- function(x){
  if (inherits(x, "sst")){
    diagels <- isondiag_vech(x[1, ])
    newx <- x / rowSums(x[, diagels])
  } else {
    tr <- sum(diag(x))
    newx <- x/tr
  }
  return(newx)
}
#' @export
normalize_trace <- normalise_trace

cov_evals_ft <- function(ms, H = NULL, evecs = NULL, av = NULL){
  if (is.null(H)){H <- helmertsub(dimfromvech(ms[1, ]))}
  H %*% cov_evals_est(ms, evecs = evecs, av = av) %*% t(H)
}

descendingordererror <- function(d0){
  # good help on withRestarts and related here: http://adv-r.had.co.nz/beyond-exception-handling.html
  withRestarts(
    stop(structure(
      class = c("est_evals_not_descending", "error", "condition"),
      list(message = paste("Estimated common eigenvalues are not in descending order:", paste(d0, collapse = " ")),
           call = sys.call(-1))
    )),
    ignore = function() d0,
    sort = function() sort(d0, decreasing = TRUE, na.last = TRUE),
    use_NA = function() NA * d0,
    use_value = function(xx) xx
  )
}
