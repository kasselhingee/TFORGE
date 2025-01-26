# Utility functions for Schwartzman 2010 test
# vecd (which is not quite vech and has a different ordering)
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
vech2vecd_mat <- function(len){
  permmat <- matrix(0, len, len)
  isdiag <- isondiag_vech(rep(1, len))
  permmat[cbind(1:sum(isdiag), which(isdiag))] <- 1
  permmat[cbind(sum(isdiag) + (1:sum(!isdiag)), which(!isdiag))] <- 1
  scalemat <- diag(c(rep(1, sum(isdiag)), rep(sqrt(2), sum(!isdiag))))
  return(scalemat %*% permmat)
}
vech2vecd <- function(vec){
  drop(vech2vecd_mat(length(vec)) %*% vec)
}

# Special sample covariance that uses Schwartzman's vecd, but with merr in vech form
S_mcovar <- function(merr){
  tmpmat <- vech2vecd_mat(ncol(merr))
  tmpmat %*% mcovar(merr) %*% t(tmpmat)
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
# These are defined in the final equation of Section 2.4 of Schwartzman 2010
#x1 and x2 are the two samples
#M1 and M2 are the means
#C1 and C2 are the covariances (Schwartzman style)
#n1 and n2 are the sample sizes
S_anv <- function(n1, n2, M1, M2, C1, C2){
  # cov1 <- S_mcovar(merr(x1))
  # cov2 <- S_mcovar(merr(x2))
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
    v = v,
    SigmaOmega_evals_av = tr_so / nrow(Sigma),
    SigmaOmega_evals2_av = tr_soso / nrow(Sigma)
  ))
}

#' @title Two Sample Test of Equal Eigenvalues Using GOE Approximation
#' @description Applies the equal-eigenvalue hypothesis test between two samples by \insertCite{schwartzman2010gr;textual}{TFORGE}.
#' The null hypothesis is that the population means of each sample have the same eigenvalues, regardless of eigenvectors.
#' The test uses a statistic from the situation that both populations are Gaussian Orthogonal Ensembles (Gaussian-distributed independent elements with the variance on the off diagonal elements half that of the diagonal elements).
#' The distribution of this statistic for more general populations is approximated using a tangent space and the Welch-Satterthwaite approximation.
#' @details 
#' The test statistic is equation 11 of \insertCite{schwartzman2010gr}{TFORGE}.
#' The p value of the test is computed using the approximate distribution reached at the end of \insertCite{@Section 2.4, @schwartzman2010gr}{TFORGE}.
#' @param x1 A single sample of matrices (passed to [`as_fsm()`]) or
#' a list of two samples of matrices (passed to [`as_kfsm()`]).
#' @param x2 If `x1` is a single sample then `x2` must be the second sample. Otherwise `x2` should be `NULL`.
#' @param scalestat If `TRUE` then the statistic divided by the estimated \eqn{a}. This modified statistic has approximately the same scale regardless of the data, although the thickness of the distribution tails (related to \eqn{v}) will vary. Simulations and general boostrap theory suggest that `scalestat=TRUE` leads to better test size and power when using bootstrap.
#' @inheritParams test_unconstrained
#' @references
#' \insertAllCited{}
#' @return
#' A list:
#' + `pval` The p-value.
#' + `t` Value of the statistic \insertCite{@Eq. 11, @schwartzman2010gr}{TFORGE}.
#' + `a` Plug-in estimate of the \eqn{a} in the final equation of \insertCite{@Section 2.4, @schwartzman2010gr}{TFORGE}.
#' + `v` Plug-in estimate of the \eqn{v} in the final equation of \insertCite{@Section 2.4, @schwartzman2010gr}{TFORGE}.
#' + `var_Lambda_evals` The variance of the eigenvalues of Schwartzman's \eqn{\Lambda}{Lambda} matrix, which may relate to the quality of the Welch-Satterthwaite approximation. 
#' @export
test_unconstrained_aGOE <- function(x, x2 = NULL, B = "chisq", nullevals = "av", scalestat = TRUE){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_kfsm")){
    if (length(x) == 2){
      if (!is.null(x2)){stop("x is a list of samples, so x2 must be NULL.")}
      x2 <- x[[2]]
      x1 <- x[[1]]
    } else if (length(x > 2)){
      stop("Must be exactly two samples")
    }
  } else {
    x2 <- as_fsm(x2)
    x1 <- as_fsm(x)
  }
  n1 <- nrow(x1)
  n2 <- nrow(x2)

  if (B == "chisq"){# use Schwartzman's asymptotic calibration
    Tstat <- stat_unconstrained_aGOE(as_flat(list(x1, x2)), scale = scalestat)
    if (scalestat){
      pval <- 1-stats::pchisq(Tstat, df = attr(Tstat, "df"))
      a = 1
      df = attr(Tstat, "df")
      var_Lambda_evals = attr(Tstat, "var_Lambda_evals")
    } else {# do the scaling of stat here
      anv <- S_anv(n1, n2, attr(Tstat, "M1"), attr(Tstat, "M2"), 
            C1 = S_mcovar(merr(x1, mean = attr(Tstat, "M1"))),
            C2 = S_mcovar(merr(x2, mean = attr(Tstat, "M2"))))
      pval <- 1-stats::pchisq(Tstat / anv$a, df = anv$v)
      a = anv$a
      df = anv$v
      var_Lambda_evals = anv$SigmaOmega_evals2_av - anv$SigmaOmega_evals_av^2
    }
    out <- list(
      pval = as.vector(pval),
      t0 = as.vector(Tstat),
      a = a,
      df = df,
      var_Lambda_evals = var_Lambda_evals
    )
    class(out) <- c("TFORGE", class(out))
    return(out)
  }
  
  # use bootstrapping from the null
  L1 <- eigen_desc(mmean(x1))$values
  L2 <- eigen_desc(mmean(x2))$values
  nullevals <- switch(nullevals,#could be anything really!
                      `1` = L1,
                      `2` = L2,
                      av = (L1 + L2)/2)
  x_std <- as_flat(list(standardise_specifiedevals(x1, nullevals), 
                        standardise_specifiedevals(x2, nullevals)))
  res <- bootresampling(as_flat(list(x1, x2)), x_std, 
                        stat = stat_unconstrained_aGOE,
                        B = B,
                        scale = scalestat)
  return(res)
}

# for similarity with stat_unconstrained
# not always scaling by a because it is computationally expensive to bootstrap
stat_unconstrained_aGOE <- function(x, scale = FALSE){
  x <- as_flat(x)
  stopifnot(inherits(x, "TFORGE_kfsm"))
  stopifnot(length(x) == 2)
  n1 <- nrow(x[[1]])
  n2 <- nrow(x[[2]])
  M1 <- mmean(x[[1]])
  M2 <- mmean(x[[2]])
  L1 <- eigen_desc(M1)$values
  L2 <- eigen_desc(M2)$values
  Tstat <- sum((L1 - L2)^2) * n1 * n2 / (n1 + n2)
  attr(Tstat, "M1") <- M1
  attr(Tstat, "M2") <- M2
  if (scale){#now for approximate scale and degrees of freedom
    anv <- S_anv(n1, n2, M1, M2, 
                 C1 = S_mcovar(merr(x[[1]], mean = M1)),
                 C2 = S_mcovar(merr(x[[2]], mean = M2)))
    Tstat <- Tstat / anv$a
    attr(Tstat, "a") <- 1
    attr(Tstat, "df") <- anv$v
    attr(Tstat, "var_Lambda_evals") <- anv$SigmaOmega_evals2_av - anv$SigmaOmega_evals_av^2 #possibly relevant to Welche-Satterthwaite approximation quality - but didn't see much of an association
  }
  return(Tstat)
}

# @title Schwartzman's theoretical Tstar statistic for eigenvalues (eq 16)
# @details The test statistic computed is (eq 16 Schwartzman 2010), using given population means M1 and M2.
# @param M1 population expectation for first population
# @param M2 population expectation for second population
statstar_schwartzman_eval <- function(x1, x2, M1, M2){
  x1 <- as_fsm(x1)
  x2 <- as_fsm(x2)
  Z1 <- mmean(x1) - M1
  Z2 <- mmean(x2) - M2
  U1 <- eigen_desc(M1)$vectors
  U2 <- eigen_desc(M2)$vectors
  Omega <- Omega_eval(nrow(x1), nrow(x2), U1, U2)
  Tstatstar <- t(c(vecd(Z1), vecd(Z2))) %*% Omega %*% c(vecd(Z1), vecd(Z2))
  return(drop(Tstatstar))
}

# simulate matrices using rmvnorm where conversion to and from vector is via vecd
rsymm_Schwartzman <- function(n, mean, sigma = diag(length(vecd(mean)))){
  stopifnot(isSymmetric(mean))
  tmp <- mvtnorm::rmvnorm(n, mean = vecd(mean), sigma = sigma)
  return(apply(tmp, MARGIN = 1, invvecd, simplify = FALSE))
}

