testthat("stat_commoneigenvals() doesn't reject for simulation of single sample from null", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  t0 <- stat_commoneigenvals(c(3, 2, 1), Ysample)
  Ystdsample <- standardise_commoneigenvals(c(3,2,1), Ysample)
  B <- 10
  nullt <- replicate(B,
    stat_commoneigenvals(c(3, 2, 1), sample(Ystdsample, replace = TRUE)))
  pval <- mean(nullt > t0)
  expect_gt(pval, 0.2)
})

testthat("stat_commoneigenvals() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  Ystdsample <- standardise_commoneigenvals(c(3,2,1), Ysample)
  expect_equal(stat_commoneigenvals(c(3,2,1), Ystdsample), 0)
})


