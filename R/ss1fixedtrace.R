stat_ss1fixedtrace <- function(x, evals = NULL){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  stopifnot(all(dim(x[[1]][[1]]) == c(3,3)))
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  
  # some building blocks
  A0 <- matrix(c( 0, 1,-1,
                 -1, 0, 1,
                  1,-1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  av <- lapply(x, mmean)
  ess <- lapply(av, eigen)
  naveval <- lapply(ess, function(a){a$values/sqrt(sum(a$values^2))})
  
  # The Omegas
  Omegas <- mapply(function(ms, av, ess){
    covar_unconstrained <- cov_evals(ms = ms, evecs = ess$vectors, av = av)
    projmat <- (diag(1, nrow(av)) - (ess$values %*% t(ess$values)/sum(ess$values^2)))/sqrt(sum(ess$values^2))
    projmat %*% covar_unconstrained %*% projmat
    },
    ms = x,
    av = av,
    ess = ess,
    SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
    #estimate according to equation after (41)
    mats <- mapply(function(aveval, Omega){
      t(A0) %*% aveval %*% t(aveval) %*% A0 / drop(t(aveval) %*% t(A0) %*% Omega %*% A0 %*% aveval)
      },
      aveval = naveval,
      Omega = Omegas,
      SIMPLIFY = FALSE)
    mat <- purrr::reduce(mats, `+`)
    eigenmat <- eigen(mat)
    idx <- max(which(eigenmat$values > sqrt(.Machine$double.eps))) #zero eigenvalue corresponds to the 1,1,1 eigenvector
    d0 <- eigenmat$vectors[, idx]
    # make d0 DOT evalsav have as much positive sign as possible
    dotprds <- vapply(naveval, function(v){v %*% d0}, FUN.VALUE = 0.1)
    avsign <- mean(sign(dotprds))
    if (avsign < 0){
      d0 <- -1 * d0
    }
    stopifnot(all(order(d0, decreasing = TRUE) == 1:length(d0)))
    warning("should estimated d0 be descending?")
  } else {
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
    if (!isTRUE(all.equal(sum(d0), sum(diag(x[[1]][[1]]))))){stop("Provided evals do not sum to trace of observations.")}
  }
  
  # now the statistic (40) for each sample:
  persamplestat <- mapply(function(d3, Omega, n){
    n * (t(d0) %*% t(A0) %*% d3)^2 / (t(d0) %*% t(A0) %*% Omega %*% A0 %*% d0)
  },
  d3 = naveval,
  Omega = Omegas,
  n = lapply(x, length),
  SIMPLIFY = FALSE
  )
  stat <- drop(purrr::reduce(persamplestat, `+`))
  attr(stat, "null_evals") <- drop(d0)
  return(stat)
}