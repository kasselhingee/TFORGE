test_that("commoneigenvals doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- replicate(5, lapply(rsymm(50, 3), `+`, diag(c(3,2,1))), simplify = FALSE) 

  t0info <- stat_commonevals_ksample(Ysamples)
  expect_equal(t0info$esteval, c(3, 2, 1), tolerance = 1E-2) 
  #est_via_av <- eigen(mmean(do.call(c, Ysamples)))$values

  Ystdsamples <- lapply(Ysamples, standardise_specifiedevals, t0info$esteval)
  tmpres <- stat_commonevals_ksample(Ystdsamples)
  expect_equal(tmpres$stat, 0)
  expect_equal(tmpres$esteval, t0info$esteval)

  B <- 10
  nullt <- replicate(B, { stat_commonevals_ksample(lapply(Ystdsamples, sample, replace = TRUE))$stat })
  pval <- mean(nullt > t0info$stat)
  expect_gt(pval, 0.2)
})
