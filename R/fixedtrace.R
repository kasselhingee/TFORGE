# functions for testing eigenvalues when the trace of the matrices are fixed

#' @title Check if the supplied sample(s) have fixed trace
#' @description Compares the trace of all the supplied matrices to check that they are equal.
#' @param x A sample or multiple samples of matrices suitable for [`as_flat()`].
#' @param tolerance Tolerance on the relative difference, passed to [`all.equal()`]
#' @export
#' @return `TRUE` or `FALSE`
has_fixedtrace <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_kfsm")){x <- do.call(rbind, x)}
  diagels <- isondiag_vech(x[1, ])
  traces <- rowSums(x[, diagels])
  tracerange <- range(traces)
  isTRUE(all.equal(tracerange[1], tracerange[2], tolerance = tolerance))
}

 
#' @title Test for eigenvalues when trace is fixed
#' @description 
#' For a single sample of symmetric matrices with fixed trace, test eigenvalues of the population mean.
#' For multiple samples of symmetric matrices with fixed trace, test for equality of the eigenvalues of the population means.
#' The test statistic is calculated by `stat_fixedtrace()`.
#' @details
#' The fixed trace constraint forces the set of eigenvalues to lie in a plane.
#' The test statistic accounts for this constraint by using an orthonormal basis in the plane.
#' Bootstrap resampling is from an empirical distribution that satisfies the null hypothesis; for this test we use empirical likelihood \insertCite{owen:2013}{TFORGE} to find non-uniform sampling weights for each matrix in the original sample.
#'
#' Eigenvalues must be distinct.
#' # Hypotheses
#' For a single sample the null hypothesis is that the population mean has eigenvalues of `evals`; the alternative hypothesis is that the eigenvalues are not equal to `evals`.
#' For multiple samples, `evals` must be omitted and the null hypothesis is that the population means of each sample have the same eigenvalues.
#' @references \insertAllCited{}
#' @inherit test_unconstrained return
# @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as_flat()`] for required structure.
# @param evals If supplied the eigenvalues of the null hypothesis. When supplied `evals` must sum to the trace of the matrices. For the multisample statistic this should be `NULL` and the null evals estimated by the function.
#' @param maxit The maximum number of Newton steps allowed in empirical likelihood optimisation \insertCite{owen:2013}{TFORGE}.
#' @inheritParams test_unconstrained
#' @export
test_fixedtrace <- function(x, evals = NULL, B, maxit = 25){
  x <- as_flat(x)
  stopifnot(has_fixedtrace(x))
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  if (!is.null(evals)){
    if (all(abs(evals - evals[1]) < sqrt(.Machine$double.eps))){
      warning("Supplied evals are equal and this test is not designed for non-distinct eigenvalues. See test_multiplicity() instead.")
    }
  }

  if (has_ss1(x)){warning("All tensors have a sum of squared eigenvalues of 1. Consider using test_ss1fixedtrace().")}
  
  if (B == "chisq"){
    df <- (dim_fsm_kfsm(x) - 1) * (length(x) - is.null(evals))
    return(chisq_calib(x, stat_fixedtrace, df = df, evals = evals))
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
    class(out) <- c("TFORGE", class(out))
    return(out)
  }

  res <- bootresampling(x, wts, 
                        stat = stat_fixedtrace,
                        B = B,
                        evals = evals)
  return(res)
}

#' @rdname test_fixedtrace
#' @details 
#' The test statistic is calculated by `stat_fixedtrace()`. 
#' @export
stat_fixedtrace <- function(x, evals = NULL){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_fsm")){mss <- as_flat(list(x))}
  else {mss <- x}
  if (is.null(evals) && (length(mss) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(mss) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}

  H <- helmertsub(dim_fsm_kfsm(mss))
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
 

#' @title Project diagonal elements to have trace of zero
#' @description
#' Projects the diagonal elements of symmetric matrices onto the plane through the origin and orthogonal to the vector (1,1,1,....,1).
#' The trace of the symmetric matrices is then zero.
#' @inheritParams test_multiplicity
#' @export
project_trace <- function(x){
  diagels <- isondiag_vech(x[1, ])
  H <- helmert(sum(diagels))
  projmat <- t(H) %*% diag(c(0,rep(1, sum(diagels)-1))) %*% H
  diags <- x[, diagels, drop = FALSE]
  x[, diagels] <- diags %*% t(projmat)
  return(x)
}
projtrace_matrix <- function(m){ #project to have trace 0
  diags <- diag(m)
  ones <- rep(1, length(diags))/sqrt(length(diags))
  newdiag <- diags - drop(diags %*% ones) * ones
  diag(m) <- newdiag
  return(m)
}

#' @title Scale symmetric matrices to have trace of one
#' @description Scales symmetric matrices by their trace, so that resulting matrices have a trace of one.
#' @inheritParams test_multiplicity
#' @details
#' The method will create `Inf` values for tensors that have a trace of zero.
#' @export
normalise_trace <- function(x){
  if (inherits(x, "TFORGE_fsm")){
    diagels <- isondiag_vech(x[1, ])
    newx <- x / rowSums(x[, diagels])
  } else {
    tr <- sum(diag(x))
    newx <- x/tr
  }
  return(newx)
}
#' @rdname normalise_trace
#' @export
normalize_trace <- normalise_trace

cov_evals_ft <- function(ms, H = NULL, evecs = NULL, av = NULL){
  if (is.null(H)){H <- helmertsub(dim_fsm_kfsm(ms))}
  H %*% cov_evals_est(ms, evecs = evecs, av = av) %*% t(H)
}

