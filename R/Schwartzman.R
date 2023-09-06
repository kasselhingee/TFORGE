# Schwartzman vecd (which is not quite vech and has a different ordering)
vecd <- function(m){
  c(diag(m), sqrt(2) * m[lower.tri(m, diag = FALSE)])
}
invvecd <- function(x){
  n <- (-1 + sqrt(8*length(x) + 1))/2
  m <- matrix(NA, nrow = n, ncol = n)
  diag(m) <- x[1:n]
  m[lower.tri(m)] <- x[(n+1):length(x)] / sqrt(2)
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

# Special sample covariance that uses Schwartzman's vecd
S_mcovar <- function(merr){
  if(is.array(merr)){if (length(dim(merr))==3){
    merr <- lapply(1:dim(merr)[3], function(i) merr[,,i])
  }}
  mcovarls <- lapply(merr, function(m){
    x <- vecd(m)
    x %*% t(x)
  })
  n <- length(merr)
  mmean(mcovarls) * n / (n-1)
}

# i is the index that correponds to 1, p is the number of rows/columns
Eii <- function(i, p){
  out <- matrix(0, p, p)
  out[i, i] <- 1
  return(out)
}

#make a block diagonal matrix from A and B
blockdiag <- function(A, B){
  out <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  rA <- nrow(A)
  cA <- ncol(A)
  out[1:rA, 1:cA] <- A
  out[rA + (1:nrow(B)), cA + (1:ncol(B))] <- B
  return(out)
}

#Schwartzman's wii from eqn 16
wii <- function(i, U1, U2){
  Emat <- Eii(i, nrow(U1))
  c(
    vecd(U1 %*% Emat %*% t(U1)),
    -vecd(U2 %*% Emat %*% t(U2))
  )
}

# Schwartzman's eigenvalue Omega(M) (eq 16)
# @param n1 Sample size of 1
# @param n2 Sample size of 2
# @param U1 Eigenvectors for mean of 1
# @param U2 Eigenvectors for mean of 2
Omega_eval <- function(n1, n2, U1, U2){
  vals <- lapply(1:nrow(U1), function(i){
    w <- wii(i, U1, U2)
    w %*% t(w)
  })
  msum(vals) * n1 * n2 / (n1 + n2)
}

#compute estimates Schwartzman's a and v for the test statistic distribution
#ms1 and ms2 are the two samples
#M1 and M2 are the means
#C1 and C2 are the covariances (Schwartzman style)
#n1 and n2 are the lengths
S_anv <- function(n1, n2, M1, M2, C1, C2){
  # cov1 <- S_mcovar(merr(ms1))
  # cov2 <- S_mcovar(merr(ms2))
  Sigma <- blockdiag(C1/n1, C2/n2)
  
  es1 <- eigen_desc(M1)
  if (is.complex(es1$values)){stop("Matrix M1 has complex eigenvalues.")}
  es2 <- eigen_desc(M2)
  if (is.complex(es2$values)){stop("Matrix M2 has complex eigenvalues.")}
  Omega <- Omega_eval(n1, n2, es1$vectors, es2$vectors)
  
  tr_soso <- sum(diag(Sigma %*% Omega %*% Sigma %*% Omega))
  tr_so <- sum(diag(Sigma %*% Omega))
  a <- tr_soso/tr_so
  v <- tr_so^2/tr_soso
  return(list(
    a = a,
    v = v
  ))
}

#' @title Schwartzman et. al. Two Sample Eigenvalue Test 
#' @description Computes the citeSchwartzman test statistic for equality between eigenvalues (eq 11 Schwartzman 2010) and corresponding approximate p value.
#' @details The test statistic computed is (eq 11 Schwartzman 2010).
#' The p value computed uses the approximate distribution reached at the end of section 2.4 Schwartzman: the distribution is chi-squared with first two moments approximating the first two moments of a statistic (eq 16 Schwartzman), that approximates (eq 11 Schwartzman 2010) for large samples.
#' @param ms1 Sample of matrices.
#' @param ms2 Sample of matrices.
#' @export
stat_schwartzman_eval <- function(ms1, ms2){
  n1 <- length(ms1)
  n2 <- length(ms2)
  M1 <- mmean(ms1)
  M2 <- mmean(ms2)
  L1 <- eigen_desc(M1)$values
  L2 <- eigen_desc(mmean(ms2))$values
  Tstat <- sum((L1 - L2)^2) * n1 * n2 / (n1 + n2)

  #now for the distribution
  anv <- S_anv(n1, n2, M1, M2, 
        C1 = S_mcovar(merr(ms1, mean = M1)),
        C2 = S_mcovar(merr(ms2, mean = M2)))
  pval <- 1-pchisq(Tstat / anv$a, df = anv$v)
  return(list(
    pval = pval,
    t = Tstat,
    a = anv$a,
    v = anv$v
  ))
}

# @title Schwartzman's theoretical Tstar statistic for eigenvalues (eq 16)
# @details The test statistic computed is (eq 16 Schwartzman 2010), using given population means M1 and M2.
# @param M1 population expectation for first population
# @param M2 population expectation for second population
#' @export
statstar_schwartzman_eval <- function(ms1, ms2, M1, M2){
  Z1 <- mmean(ms1) - M1
  Z2 <- mmean(ms2) - M2
  U1 <- eigen_desc(M1)$vectors
  U2 <- eigen_desc(M2)$vectors
  Omega <- Omega_eval(length(ms1), length(ms2), U1, U2)
  Tstatstar <- t(c(vecd(Z1), vecd(Z2))) %*% Omega %*% c(vecd(Z1), vecd(Z2))
  return(drop(Tstatstar))
}

# simulate matrices using rmvnorm where conversion to and from vector is via vecd
rsymm_Schwartzman <- function(n, mean, sigma = diag(length(vecd(mean)))){
  stopifnot(isSymmetric(mean))
  tmp <- mvtnorm::rmvnorm(n, mean = vecd(mean), sigma = sigma)
  return(apply(tmp, MARGIN = 1, invvecd, simplify = FALSE))
}

