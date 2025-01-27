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
