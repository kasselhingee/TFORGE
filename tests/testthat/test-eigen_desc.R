test_that("eigen_desc gives descending eigenvalues and matching vectors", {
  set.seed(345)
  m <- inv_vech(rsymm_norm(1, diag(c(1,0,-2)), sigma = 0.1*diag(1, 3*4/2))[1, ])
  ess <- eigen_desc(m)
  
  expect_equal(ess$vectors %*% diag(ess$values) %*% t(ess$vectors), m)
  expect_equal(order(ess$values, decreasing = TRUE), 1:3)
})
