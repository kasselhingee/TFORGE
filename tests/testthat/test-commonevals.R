test_that("test_commonevals() doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- replicate(5, lapply(rsymm(50, 3), `+`, diag(c(3,2,1))), simplify = FALSE) 

  res <- test_commonevals(Ysamples, 20)
  expect_gt(res$pval, 0.2)
})

test_that("test_commonevals() repeats with given .Random.seed value", {
  set.seed(34)
  Ysamples <- replicate(5, lapply(rsymm(50, 3), `+`, diag(c(3,2,1))), simplify = FALSE) 
  res1 <- test_commonevals(Ysamples, 10)
  .Random.seed <- res1$seed
  res2 <- test_commonevals(Ysamples, 10)
  expect_equal(res1, res2)
})

test_that("test_commonevals() behaves elegantly when evals can't be estimated", {
  set.seed(2134)
  mss <- list(lapply(rsymm(10, 3, 1), `+`, diag(c(3,2,1))),
              lapply(rsymm(10, 3, 1), `+`, diag(c(4,3,2))))
  expect_warning(res <- test_commonevals(mss, 100), "could not be calculated")
})
