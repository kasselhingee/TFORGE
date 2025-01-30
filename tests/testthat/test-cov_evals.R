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

