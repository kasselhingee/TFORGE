test_that("cov_evals() and cov_evals_old() match", {
  set.seed(12)
  Ysample <- as_fsm(rsymm(1E2, diag(c(3,2,1))))
  mcov <- mcovar(Ysample)
  
  # any vectors:
  set.seed(4)
  evecs <- runifortho(3)
  expect_equal(cov_evals(evecs, mcov), cov_evals_old(evecs, mcov))
  
  # eigenvectors:
  evecs <- eigen_desc(mmean(Ysample))$vectors
  expect_equal(cov_evals(evecs, mcov), cov_evals_old(evecs, mcov))
})

