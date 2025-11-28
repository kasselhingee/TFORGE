#' @title Test eigenvalue multiplicity
#' @description
#' Tests the multiplicity of the eigenvalues a population's mean.
#' The test statistic is computed by `stat_multiplicity()`.
#' The null hypothesis is that the population mean has the specified the multiplicity of eigenvalues.
#' For unconstrained symmetric matrices or symmetric matrices with fixed trace use `test_multiplicity()`. 
#' For matrices constrained to have non-negative eigenvalues use `test_multiplicity_nonnegative()`.
#' @details
#' For `test_multiplicity()`, bootstrap resampling is conducted from the null hypothesis by first translating the original sample to satisfy the null hypothesis with `translate2multiplicity()`. 
#' For `test_multiplicity_nonnegative()`, weighted bootstrapping is used (see 'Weighted Bootstrapping' below).
#'
#' On `refbasis`:
#' An estimate of each eigenspace specified by `mult` can be obtained from the eigenvectors of the sample mean.
#' The eigenvectors create an orthonormal basis of the (estimated) eigenspace,
#' however the choice of orthonormal basis for the estimated eigenspace effects the performance.
#' This choice is specified by the parameter `refbasis`.
#' Setting `refbasis = "sample"` uses the eigenvectors of the sample mean as the basis, however the resulting statistic does not appear to be pivotal.
#' Choosing the orthonormal basis independently of the data does result in a pivotal asymptotically chi-squared statistic.
#' Setting `refbasis = "random"` will do exactly this, by applying a uniformly random rotation of the relevant eigenvectors of the sample mean.
#' We recommend using `refbasis = "sample"` (which requires bootstrap calibration) because test power is much higher than `refbasis = "random"`.
#' We recommend that the number bootstrap resamples is at least 1000 if `refbasis = "sample"`.
#'
#' For 3x3 matrices, the weighted-bootstrapping method used by `test_multiplicity_nonnegative()` has poor test size for samples smaller than 20;
#' larger matrices will likely need larger samples.
#' 
#' Due to the random rotation of the eigenvectors when `refbasis = "random"`, use [`set.seed()`] if you want the answer to be repeatable.
#' @param x A sample of matrices suitable for [`as_fsm()`].
#' @param mult A vector specifying the eigenvalue multiplicity under the null hypothesis in descending order of eigenvalue size.
#' @param refbasis Select the basis of the eigenspaces. See details.
#' @inheritParams test_unconstrained
#' @inheritParams test_fixedtrace
#' @inheritSection test_fixedtrace Weighted Bootstrapping
#' @examples
#' x <- rsymm_norm(15, mean = diag(c(2, 1, 1, 0)))
#' test_multiplicity(x, mult = c(1, 2, 1))
#' @return A `TFORGE` object (see [`boot_calib()`] or [`chisq_calib()`]) including p-value of the test (slot `pval`) and the statistic for `x` (slot `t0`). 
#' @export
test_multiplicity <- function(x, mult, B = 1000, refbasis = "sample"){
  x <- as_flat(x)
  if (B == "chisq"){
    if (refbasis[[1]] == "sample"){warning("chisq calibration does not work for the statistic when using the eigenvectors of the sample mean (refbasis = 'sample')")}
    return(chisq_calib(x, stat_multiplicity, df = sum(mult) - length(mult), mult = mult, refbasis = refbasis))
  }
  ms_std <- translate2multiplicity(x, mult)
  res <- boot_calib(x, ms_std, 
    stat = stat_multiplicity,
    B = B,
    mult = mult,
    refbasis = refbasis)
  return(res)
}

#' @rdname test_multiplicity
#' @export
test_multiplicity_nonnegative <- function(x, mult, B = 1000, maxit = 25, refbasis = "sample"){
  x <- as_flat(x)
  if (B == "chisq"){
    if (refbasis == "sample"){warning("chisq calibration does not work for the statistic using eigenvalues of the sample mean")}
    return(chisq_calib(x, stat_multiplicity, df = sum(mult) - length(mult), mult = mult, refbasis = refbasis))
  }
  
  # compute corresponding weights that lead to maximum emp.lik.
  av <- mmean(x)
  nullmean <- multiplicity_nullmean(av, mult)
  wts <- elwts_fixedtrace(x, nullmean, maxit)
  
  #check the weights
  if (!wtsokay(wts)){
    out <- list(
      pval = 0,
      t0 = stat_multiplicity(x, mult, refbasis = refbasis),
      nullt = NA,
      stdx = wts,
      B = NA
    )
    class(out) <- c("TFORGE", class(out))
    return(out)
  }
  
  res <- boot_calib(x, wts, 
                        stat = stat_multiplicity,
                        B = B,
                        mult = mult,
                        refbasis = refbasis)
  return(res)
}

#' @rdname test_multiplicity
#' @param evecs For debugging only. Supply eigenvectors of population mean.
#' @export
stat_multiplicity <- function(x, mult, evecs = NULL, refbasis = "sample"){
  av <- mmean(x)
  if (sum(mult) != ncol(av)){
    stop(paste("Sum of mult = ", mult, "is not equal to ", ncol(av), collapse = " "))
  }
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  if (!is.null(evecs)){
    warning("evecs should only be supplied for debugging")
    es <- list()
    es$vectors <- evecs
  } else {
    es <- eigen_desc(av)
  }
  C0 <- mcovar(merr(x, mean = av)) # the covariance between elements of xi

  #indices of multiplicity
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  
  # Below chooses the basis for each eigenspace. 
  # Arbitrary options are 'random' or projection-like from a set of prespecified eigenvectors (e.g. cannonical vectors).
  # Non-arbitrary are 'sample' - uses the eigenvectors of the sample mean,
  # and 'mincorr' which finds eigenvectors that minimises the correlation between the resulting eigenvalue estimates
  if (is.character(refbasis) && (refbasis[[1]] == "sample")){
    es$vectors <- es$vectors
  } else if (is.character(refbasis) && (refbasis[[1]] == "mincorr")){
    es$vectors <- arbitrary_evecs(es$vectors, idxs, refbasis = diag(nrow(es$vectors))) #start with canonical vectors
    es$vectors <- mostdiagevecs(es$vectors, idxs, mcov = C0/nrow(x)) #optimise vectors to reduce correlation
  } else {
    es$vectors <- arbitrary_evecs(es$vectors, idxs, refbasis = refbasis) #project-like specified vectors to become orthogonal bases of the estimated eigenspaces
  }
  es$values <- diag(t(es$vectors) %*% av %*% es$vectors)
  
  # for each unique hypothesised eigenvalue, projects the estimated eigenvalues into the plane perpendicular to (1,1,...1)
  xi <- xiget(es$values, mult, idxs)

  covar <- xicovar(mult, idxs, es$vectors, C0/nrow(x))

  return(drop(t(xi) %*% solve_error(covar) %*% xi))
}

#this function turns the qYbarq into vectors perpendicular to 1,1,1,1,1,1...
evalsvs111 <- function(evals){
    if (length(evals) < 2){return(NULL)}
    wmat <- helmertsub(length(evals))
    return(wmat %*% evals)
}

# the random variables xi in sets per multiplicity because the weight matrix is different in each one
# @param idxs gives the locations in the full eigenvalue vector of the eigenvalues corresponding to each multiplicity
xiget <- function(evals, mult, idxs){
  xi <- lapply(1:length(mult), function(j){
    evals <- evals[ idxs[[j]] ]
    return(evalsvs111(evals))
  })
  xi <- unlist(xi)
  return(xi)
}

# compute/estimate covariance of xi, Cav is the covariance of av = C0/n
xicovar <- function(mult, idxs, evecs, Cav){
  # for each pair of distinct eigenvalues with multiplicity>1, create the matrix of covariance between xi, these are entries in the block-matrix form of V3 in my notes 
  valpairs <- expand.grid(1:length(idxs), 1:length(idxs))
  blocks <- apply(valpairs, MARGIN = 1, function(jk){
    if (mult[jk[1]] < 1.5){return(NULL)}
    if (mult[jk[2]] < 1.5){return(NULL)}
    helmertsub(mult[jk[1]]) %*% covarbetweenevals(jk[1], jk[2], idxs, evecs, Cav) %*% t(helmertsub(mult[jk[2]]))
  }, simplify = FALSE)

  #bind the blocks together appropriately
  # bind by row, columns still seperate
  cblocks <- lapply(1:length(mult), function(k){
     purrr::reduce(blocks[valpairs[, 2] == k], `rbind`)
  })
  covar <- purrr::reduce(cblocks, `cbind`)
  return(covar)
}

# for each pair of distinct eigenvalues with multiplicity>1, create the matrix of covariance between eigenvalues (I've called this matrix A_jk in my notes)
# j and k refer to index of *distinct* eigenvalues, i.e. entries of mult
# idxs is a list of vectors, each vector contains the columns of evecs that have columns corresponding to the vectors eigenvalue
# the computations here overlap a lot with the computation of cov_evals()
covarbetweenevals <- function(j, k, idxs, evecs, Cav){
  Dp <- dup(nrow(evecs))
  idxj = idxs[[j]]
  idxk = idxs[[k]]
  if ((length(idxj) == 1 ) || (length(idxk) == 1)){return(NULL)}
  vecpairs <- expand.grid(idxj, idxk)
  At <- apply(vecpairs, MARGIN = 1, function(uv){
    # qjuqkv <- evecs[, ] %*% t(evecs[, uv[2]])
    qjuqju <- evecs[, uv[1]] %*% t(evecs[, uv[1]])
    qkvqkv <- evecs[, uv[2]] %*% t(evecs[, uv[2]])
    t(vec(qjuqju)) %*% Dp %*% Cav %*% t(Dp) %*% vec(t(qkvqkv))
  }, simplify = TRUE)
  dim(At) <- c(length(idxj), length(idxk)) # converts to matrix taking advantage that the first index is filled fastest
  return(At) # this is actually the transpose of whats in my notes (where rows correspond to idxk)
}

# Schwartzman's blk() operation for eigenvalues
# Returns a modification of evals that matches multiplicity given in mult
# Assumes evals in descending order, and 
# averages in these evals blocks given by mult
# Useful for computing a hypothetical population mean that satisfies the null hypothesis
multiplicity_blk <- function(evals, mult){
  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  
  
  # get eigenvalues for each multiplicity by averaging 
  newevals <- lapply(1:length(mult), function(i){
    newval <- mean(evals[esvalstart[i] : cmult[i]])
    rep(newval, mult[i])
  })
  newevals <- unlist(newevals)
  return(newevals)
}

# From the sample average, calculate a hypothetical population mean that satisfies the given hypothetical multiplicity
multiplicity_nullmean <- function(av, mult){
  stopifnot(sum(mult) == ncol(av))
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  es <- eigen_desc(av, symmetric = TRUE)
  
  # evals according to null
  nullevals <- multiplicity_blk(es$values, mult)
  
  # add in eigenvectors of av to get full null mean
  nullmean <- es$vectors %*% diag(nullevals) %*% t(es$vectors)
  return(nullmean)
}

#' @rdname test_multiplicity
#' @export
translate2multiplicity <- function(x, mult){
  x <- as_fsm(x)
  av <- mmean(x)
  nullmean <- multiplicity_nullmean(av, mult)
  
  # make a translated sample out of nullmean and errors
  nullmean <- vech(nullmean)
  av <- vech(av)
  out <- t(t(x) - av + nullmean)
}

# For random rotations of the eigenspace basis, want random orthogonal matrices with positive determinant.
# Here I dont bother with the positive determinant requirement as the sign of the axes doesnt matter.
# Random matrices in the uniform distribution on the Stiefel manifold can be obtained as the L in X = TL decomposition of a random Normal matrix (elements iid normal), where
# T is lower triangular (and square) and L is unitary (not necessarily square).
# [P67-68 Gupta and Nagar 1999]
# This is not quite the QR decomposition in LINPACK which is X = QR such that Q is orthonormal (square matrix) and R is upper triangular.
# t(X) = QR => X = R^T Q^T.
# R^T is lower triangular ==> T
# Q^T is orthonormal ==> L.
# So when X is square the Q of the QR decomposition of t(X) is the transpose of a uniformly random matrix on the Stiefel manifold.
# @param p the nrows and columns
runifortho <- function(p){
  # rstiefel::rustiefel(p,p)
  m <- matrix(stats::rnorm(p * p),
         nrow = p,
         ncol = p)
  out <- qr.Q(qr(t(m), LAPACK = TRUE))
  return(t(out))
}

#' @noRd
#' @title randomly rotate eigenvectors for each eigenspace
#' @description Get eigenvalue estimates after arbitrary choice of the basis defining the eigenspace
# To uniformly randomly rotate them in each block, need uniform rotation matrices using runifortho()
# For each eigenspace, calls arbitrary_basis() to do the actual calculations.
arbitrary_evecs <- function(evecs, idxs, refbasis = "random"){
  #check refbasis:
  if (is.character(refbasis)){if (refbasis != "random"){stop("refbasis must be either 'random' or a matrix")}}
  else{
    stopifnot(is.matrix(refbasis))
    stopifnot(nrow(refbasis) == nrow(evecs))
    stopifnot(nrow(refbasis) == ncol(refbasis))
    if (!all(abs(refbasis %*% t(refbasis) - diag(1, nrow(refbasis))) < sqrt(.Machine$double.eps))){stop("refbasis does not appear to be orthogonal")}
  }
  
  # Do calculations
  newevecs <- lapply(idxs, function(idxforeval){
    if (length(idxforeval) == 1){return(evecs[, idxforeval, drop = FALSE])}
    arbitrary_basis(evecs[, idxforeval, drop = FALSE], refbasis = refbasis)
  })
  do.call(cbind, newevecs)
}


#' @noRd
#' @description For a subspace, uniformly randomly rotate the basis vectors within the subspace, or create a new basis using a reference
#' @param subspace A matrix of column vectors.
#' @param refbasis A matrix of column vectors or "random"
arbitrary_basis <- function(subspace, refbasis = "random"){
  if (is.character(refbasis)){
    if (refbasis == "random"){
    rot <- runifortho(ncol(subspace))
    return(t(rot %*% t(subspace)))
    }
  } else {
    return(project_basis(subspace, refbasis = refbasis))
  }
}

#' @noRd
#' @title Create an orthonormal basis to a lower dimensional space using a reference basis for the full space
#' @param subspace Column vectors forming an orthonormal basis of the subspace
#' @param refbasis Column vectors forming an orthonormal basis of the full space
project_basis <- function(subspace, refbasis = diag(nrow = nrow(subspace))) {
  # projection matrix is property of subspace, not the vectors used to represent it, so projmat doesn't depend on arbitrary estimation of subspace basis vectors
  projmat <- subspace %*% t(subspace)
  p <- nrow(subspace)
  # below scales refbasis axes so that alignment occurs
  # suppose x is expressed wrt ref basis
  # 1) scale values of x to make the problem like finding an ellipse
  # 2) convert x to cannonical basis by refbasis %*%
  # 3) project x to subspace
  alignedsvd <- svd(projmat %*% refbasis %*% diag(p:1), nu = ncol(subspace), nv = ncol(subspace))
  
  newbasis <- alignedsvd$u
  
  # check new basis
  stopifnot(all(abs(newbasis %*% t(newbasis) - projmat) < sqrt(.Machine$double.eps)))
  return(newbasis)
}

# find evecs that minimise the mixing of variance
mostdiagevecs <- function(evecs, idxs, mcov){
  # Do calculations
  newevecs <- lapply(idxs, function(idxforeval){
    if (length(idxforeval) == 1){return(evecs[, idxforeval, drop = FALSE])}
    d <- length(idxforeval)
    bestpar <- stats::optim(par = rep(0, (d-1)*d/2),
                  fn = ssqoffdiagonal,
                  baseevecs = evecs[, idxforeval],
                  mcov = mcov,
                  method = if (d==2){"Brent"}else{"Nelder-Mead"}, #BFGS is very similar time wise in a simulation
                  lower = if (d==2){-100}else{-Inf},
                  upper = if (d==2){100}else{+Inf},
                  control = list(maxit = 100, reltol = 1E-2)
    )
    bestrot <- cayleyTransform(vecskewsym_inverse(bestpar$par))
    bestevecs <- evecs[, idxforeval] %*% bestrot
    return(bestevecs)
  })
  do.call(cbind, newevecs)
}

vecskewsym <- function(A){A[lower.tri(A)]}
vecskewsym_inverse <- function(vec){
  p <- (1 + sqrt(8*length(vec) + 1))/2
  A <- matrix(0, p, p)
  A[lower.tri(A)] <- vec
  A[upper.tri(A)] <- -t(A)[upper.tri(A)]
  return(A)
}

#baseevecs specify the eigenspace
#vec chooses rotations within the space
#covariance of main matrix
ssqoffdiagonal <- function(vec, baseevecs, mcov){
  A <- vecskewsym_inverse(vec)
  rotmat <- cayleyTransform(A)
  trialevecs <- baseevecs %*% rotmat
  evalcov <- cov_evals(trialevecs, mcov) # cov_evals2 seems very slightly slower in the current implementation
  sum(stats::cov2cor(evalcov)[lower.tri(evalcov)]^2)
}

#A must be skew symmetric
cayleyTransform <- function(A){
  solve(diag(nrow(A)) - A) %*% (diag(nrow(A)) + A)
}
# M must be orthogonal with determinant of +1 (or zero?)
inverseCayleyTransform <- function(M){
  (M - diag(nrow(M))) %*% solve(M + diag(nrow(M)))
}
