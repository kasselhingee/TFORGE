#' @title Simulate from a matrix-t distribution
#' @details Simulates two independent matrices from matrix-Gaussian distributions and then combines them. The first matrix, `Y` say, is simulated with a covariance of `solve(Sigma) %x% diag(n + p -1)`, so that `Y %*% t(Y)` follows (according to Gupta and Nager 1999, theorem 3.2.2) a Wishart distribution with `n+p-1` degrees of freedom and scale `solve(Sigma)`.
#' The second matrix, `X` say, is simulated with a covariance of `diag(p) %x% Omega`.
#' Theorem 4.3.1 (Gupta and Nager, 1999) then says that `t(solve(Y)) %*% X + M` has as a matrix-T distribution with `n` degrees of freedom, mean `M`, and spread `Sigma` and `Omega`.
#' Returns matrices with `p` rows and `m` columns.
#' Only works when `n > p`, otherwise will need to use Cholesky factorisation

#' @param N Number of samples to draw.
#' @param n degrees of freedom, must be larger than number of rows of final matrix
#' @param M Mean matrix. Dimension p x m.
#' @param Sigma A p x p spread matrix, positive definite
#' @param Omega An m x m spread matrix, positive definite
#' @examples
#' p <- 3
#' m <- 4
#' n <- 10
#' Sigma <- diag(1, p)
#' Omega <- diag(1, m)
#' M <- matrix(0, nrow = p, ncol = m)
#' rmatrixt(2, n, M, Sigma, Omega)
#' @export
rmatrixt <- function(N, n, M, Sigma, Omega){
  p <- nrow(Sigma)
  # first simulate Y s.t. S = Y%*%t(Y) = W(n + p - 1, inv(Sigma))
  # which means Y ~ N(0, inv(Sigma) %kronecker% I_{n + p - 1}) (theorem 3.2.2)
  Ys <- mvtnorm::rmvnorm(N, mean = rep(0, p*(n+p-1)), sigma = solve(Sigma) %x% diag(n + p -1)) |>
    apply(MARGIN = 1, invvec, nrow = p, simplify = FALSE)
  Ss <- lapply(Ys, function(Y){Y %*% t(Y)})
  
  # now simulate X ~ N(0, I_p %x% Omega) (also from theorem 4.3.1)
  Xs <- mvtnorm::rmvnorm(N, mean = rep(0, p*m), sigma = diag(p) %x% Omega) |>
    apply(MARGIN = 1, invvec, nrow = p, simplify = FALSE)
  matTs <- mapply(function(S, X, M){t(solve(S)) %*% X + M}, S = Ss, X = Xs, MoreArgs = list(M = M), SIMPLIFY = FALSE)
  return(matTs)
}
