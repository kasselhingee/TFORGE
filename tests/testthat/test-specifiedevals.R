test_that("stat_commoneigenvals() doesn't reject for simulation of single sample from null", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  Ystdsample <- standardise_commoneigenvals(c(3,2,1), Ysample)
  res <- singlesampletest(Ysample, Ystdsample, 
    stat = function(Y){stat_commoneigenvals(c(3, 2, 1), Y)},
    B =  100)
  expect_gt(res$pval, 0.2)
})

test_that("stat_commoneigenvals() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  Ystdsample <- standardise_commoneigenvals(c(3,2,1), Ysample)
  expect_equal(stat_commoneigenvals(c(3,2,1), Ystdsample), 0)
})

