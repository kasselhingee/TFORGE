
test_that("a simple null situation is correctly not rejected", {
  set.seed(13)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) {m/sum(diag(m))})
    Y
  }, simplify = FALSE)
  
  stat_ss_fixedtrace(Ysamples[[1]], evals = c(3,2,1)/6)
  stat_ss_fixedtrace(Ysamples[[1]], evals = c(3,2,1)/6, evecs = diag(1, nrow = 3))
})


test_that("hasfixedtrace() gives TRUE or FALSE values", {
  set.seed(13)
  const <- 1.45
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y <- lapply(Y, function(m) {const * m/sum(diag(m))})
    Y
    }, simplify = FALSE)

  expect_true(hasfixedtrace(Ysamples))
  expect_true(hasfixedtrace(Ysamples[[1]]))
  
  set.seed(134)
  Ysamples <- replicate(5, {
    Y <- rsymm_norm(50, diag(c(3,2,1)))
    Y
  }, simplify = FALSE)
  
  expect_false(hasfixedtrace(Ysamples))
  expect_false(hasfixedtrace(Ysamples[[1]]))
})
