
# Schwartzmann vecd (which is not quite vech and has a different ordering)
vecd <- function(m){
  c(diag(m), sqrt(2) * m[lower.tri(m, diag = FALSE)])
}

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

#Schwartzmann's wii from eqn 16
wii <- function(i, U1, U2){
  Emat <- Eii(i, nrow(U1))
  c(
    vecd(U1 %*% Emat %*% t(U1)),
    -vecd(U2 %*% Emat %*% t(U2))
  )
}

#Schwartzmann's eigenvalue Omega(M) (eq 16)
#' @param n1 Sample size of 1
#' @param n2 Sample size of 2
#' @param U1 Eigenvectors for mean of 1
#' @param U2 Eigenvectors for mean of 2
Omega_eval <- function(n1, n2, U1, U2){
  vals <- lapply(1:nrow(U1), function(i){
    w <- wii(i, U1, U2)
    w %*% t(w)
  })
  msum(vals) * n1 * n2 / (n1 + n2)
}

#compute estimates Schwartzmann's a and v for the test statistic distribution
#ms1 and ms2 are the two samples
#M1 and M2 are the means
#C1 and C2 are the covariances (Schwartzmann style)
#n1 and n2 are the lengths
S_anv <- function(n1, n2, M1, M2, C1, C2){
  # cov1 <- S_mcovar(merr(ms1))
  # cov2 <- S_mcovar(merr(ms2))
  Sigma <- blockdiag(C1, C2)
  
  es1 <- eigen(M1)
  es2 <- eigen(M2)
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

stat_schwartzmann_eval <- function(ms1, ms2){
  n1 <- length(ms1)
  n2 <- length(ms2)
  M1 <- mmean(ms1)
  M2 <- mmean(ms2)
  L1 <- eigen(M1)$values
  L2 <- eigen(mmean(ms2))$values
  Tstat <- sum((L1 - L2)^2) * n1 * n2 / (n1 + n2)

  #now for the distribution
  anv <- S_anv(n1, n2, M1, M2, 
        C1 = S_mcovar(merr(ms1, mean = M1)),
        C2 = S_mcovar(merr(ms2, mean = M2)))
  pval <- pchisq(Tstat / anv$a, df = anv$v)
  return(list(
    pval = pval,
    t = tstat,
    a = anv$a,
    v = anv$v
  ))
}


