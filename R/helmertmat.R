#' @title Create the Helmert matrix
#' @param p Dimension
helmert <- function(p){
  stopifnot(p>1)
  m <- diag(c(1, -(1:(p-1))))
  m[1, ] <- 1
  m[lower.tri(m, diag = FALSE)] <- 1

  #make rows of size 1
  m <- m/sqrt(rowSums(m^2))
  return(m)
}

helmertsub <- function(p){
  helmert(p)[-1, , drop = FALSE]
}

