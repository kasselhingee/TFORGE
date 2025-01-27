test_that("cov_evals_inside() matches cov_evals_inside_cpp()",{
  set.seed(12)
  Ysample <- as_fsm(rsymm(1E2, diag(c(3,2,1))))
  mcov <- mcovar(Ysample)
  evecs <- eigen_desc(mmean(Ysample))$vectors
  dupmat <- dup(nrow(evecs))
  expect_equal(cov_evals_inside(evecs[, 1], evecs[, 2], dupmat, mcov), 
    cov_evals_inside_cpp(evecs[, 1], evecs[, 2], dupmat, mcov))
})

test_that("cov_evals() and cov_evals2() match", {
  set.seed(12)
  Ysample <- as_fsm(rsymm(1E2, diag(c(3,2,1))))
  mcov <- mcovar(Ysample)
  
  # any vectors:
  set.seed(4)
  evecs <- runifortho(3)
  expect_equal(cov_evals(evecs, mcov), cov_evals2(evecs, mcov))
  
  # eigenvectors:
  evecs <- eigen_desc(mmean(Ysample))$vectors
  expect_equal(cov_evals(evecs, mcov), cov_evals2(evecs, mcov))
})

covS <- matrix(0, 6, 6)
diag(covS) <- 1
covS[c(1,4,6), c(1,4,6)] <- 
  matrix(c(1, 0, 0,
           0, 1, 0.8,
           0, 0.8, 1), ncol = 3, nrow = 3, byrow = 3)
evalcov <- cov_evals(diag(3), covS)
diagonalizer <- eigen_desc(evalcov)$vectors
evecxevec <- apply(diag(3), 2, FUN = function(v){kronecker(v, v)})
evecxevec_new <- evecxevec %*% diagonalizer
  
svd(evecxevec_new)$u #these eigenvectors are looking hard to get from a kronecker
idx <- c(1,5,9)
a <- sqrt(svd(evecxevec_new)$u[9,3])
evec2 <- cbind(c(1,0,0), c(0,a,a), c(0,a,-a))
apply(evec2, 2, FUN = function(v){kronecker(v, v)})
# maybe evec_new isnt possible
# Is there something in the literature on getting maximally diagonal?
# Ledoit and Wolf (2004)?

ssqoffdiagonal <- function(vec, mcov){
  A <- matrix(0, 3, 3)
  A[lower.tri(A)] <- vec
  A[upper.tri(A)] <- -t(A)[upper.tri(A)]
  trialevecs <- sphm:::cayleyTransform(A)
  evalcov <- cov_evals2(trialevecs, mcov)
  sum(evalcov[lower.tri(A)]^2)
}
best <- optim(par = sphm:::inverseCayleyTransform(diag(1, 3))[lower.tri(diag(1, 3))],
      fn = ssqoffdiagonal,
      mcov = covS
      )
A <- matrix(0, 3, 3)
A[lower.tri(A)] <- best$par
A[upper.tri(A)] <- -t(A)[upper.tri(A)]
bestevecs <- sphm:::cayleyTransform(A)
cov_evals(bestevecs, covS)
