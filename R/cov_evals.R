#' Compute the Covariance of Eigenvalues
#' @description
#' For a random symmetric matrix `Y`, calculates the covariance of the eigenvalues of `Y` using the covariance of the elements of `Y` and the eigenvectors of the mean of `Y`.
#' @param evecs Matrix with columns that are eigenvectors of the mean of `Y`.
#' @param mcov Covariance of `vech(Y)`, where `Y` is the random matrix.
#' @details
#' For any two columns \eqn{a} and \eqn{b} of `evecs`, computes the covariance
#' \deqn{
#' \textrm{Cov}(a^\top Y a, b^\top Y b) = ( a \otimes a)^\top \mathbb{D} C_0 \mathbb{D}^\top (b \otimes b),
#' }
#' where \eqn{a} and \eqn{b} are the columns of `evecs` and \eqn{C_0}=`mcov` is the covariance of vech\eqn{(Y)}. \eqn{\mathbb{D}} and \eqn{\otimes} is the duplication matrix and kronecker product respectively.
#' 
#' The returned matrix has rows and columns that are in the same order as the columns of `evecs`.
#'
#' When the eigenvalues are distinct, then passing estimated eigenvectors to `cov_evals()` yields an estimate of the asymptotic covariance of the eigenvalues.
#' 
#' See Supplement B.2 for more information and derivation.
#' Uses `RcppArmadillo` for a slight increase in speed.
#' @useDynLib TFORGE, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @return A symmetric matrix with same number of columns as `evecs`.
cov_evals_old <- function(evecs, mcov){
  # the V1 / V2 matrix for a single sample depending on whether evecs estimated or supplied
  indx <- rbind(t(utils::combn(1:ncol(evecs), 2)), #this avoids repeating elements that are symmetric
                cbind(1:ncol(evecs), 1:ncol(evecs)))
  dupmat <- dup(nrow(evecs)) # is sparse so could be even faster
  vals <- mapply(cov_evals_inside_cpp, 
                 lapply(indx[,1], function(i) {evecs[, i]}),
                 lapply(indx[,2], function(i) {evecs[, i]}),
                 MoreArgs = list(dupmat = dupmat, mcov = mcov))
  V <- matrix(NA, nrow = ncol(evecs), ncol = ncol(evecs))
  V[as.matrix(indx)] <- vals
  V[lower.tri(V)] <- t(V)[lower.tri(V)]
  return(V)
}

# easier to interpret formula that gives the same results as _old:
#' @export
cov_evals <- function(evecs, mcov){
  dupmat <- dup(nrow(evecs)) # is sparse so could be even faster
  evecxevec <- apply(evecs, 2, FUN = function(v){kronecker(v, v)}) #1, 5, 9
  t(evecxevec) %*% dupmat %*% mcov %*% t(dupmat) %*% evecxevec
}


# @param x A single sample (list of matrices)
# @param av To save computation, `av` can be passed if it is already computed.
# @details If `evecs` is not provided then the eigenvectors \eqn{q_{0i}} are replaced with the eigenvectors of the average of `x`.
# @return An estimated covariance matrix for the eigenvalues of `x`.
cov_evals_est <- function(x, evecs = NULL, av = NULL){
  x <- as_flat(x)
  stopifnot(inherits(x, "TFORGE_fsm"))
  if (is.null(av)){av <- mmean(x)}
  if (is.null(evecs)){
    evecs <- eigen_desc(av, symmetric = TRUE)$vectors
  }
  mcov <- mcovar(merr(x, mean = av))
  if (all(mcov == 0)){
    stop(structure(
      class = c("zerocovariance", "error", "condition"),
      list(message = "covariance is zero",
           call =  sys.call(0))
    ))
  }

  cov_evals(evecs, mcov)
}

cov_evals_inside <- function(vecj, veck, dupmat, mcov){
  tmp <- vecj %*% t(veck)
  sum(diag(t(dupmat) %*% kronecker(tmp, tmp) %*% dupmat %*% mcov))
}


