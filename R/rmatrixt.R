#' @title Simulate from a matrix-t distribution
#' @details Simulates two independent matrices from Gaussian distributions and then combines.
#' Following notation from Gupter and Naper.
#' Returns matrices with `p` rows and `m` columns.
#' Only works when `n > p`, otherwise will need to use Cholesky factorisation

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
rmatrixt <- function(n, M, Sigma, Omega){
  p <- nrow(Sigma)
  # first simulate Y s.t. S = Y%*%t(Y) = W(n + p - 1, inv(Sigma))
  # which means Y ~ N(0, inv(Sigma) %kronecker% I_{n + p - 1}) (theorem 3.2.2)
  Y <- mvtnorm::rmvnorm(1, mean = rep(0, p*(n+p-1)), sigma = solve(Sigma) %x% diag(n + p -1)) |>
    drop() |> invvec(nrow = p)
  S <- Y %*% t(Y)
  
  # now simulate X ~ N(0, I_p %x% Omega) (also from theorem 4.3.1)
  X <- mvtnorm::rmvnorm(1, mean = rep(0, p*m), sigma = diag(p) %x% Omega) |>
    drop() |> invvec(nrow = p)
  matT <- t(solve(S)) %*% X + M
  return(matT)
}
