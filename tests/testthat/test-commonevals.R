test_that("test_commonevals() doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)

  res <- test_commonevals(Ysamples, 20)
  expect_gt(res$pval, 0.2)
})

test_that("test_commonevals() repeats with given .Random.seed value", {
  set.seed(34)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE) 
  res1 <- test_commonevals(Ysamples, 10)
  .Random.seed <- res1$seed
  res2 <- test_commonevals(Ysamples, 10)
  expect_equal(res1, res2)
})
