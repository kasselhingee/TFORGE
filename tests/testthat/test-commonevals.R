test_that("commoneigenvals doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- replicate(5, lapply(rsymm(50, 3), `+`, diag(c(3,2,1))), simplify = FALSE) 

  res <- test_commonevals(Ysamples, 100)
  expect_gt(res$pval, 0.2)
})
