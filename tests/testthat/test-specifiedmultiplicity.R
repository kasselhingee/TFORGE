test_that("stat_specifiedmultiplicity() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  av <- mmean(Ysample)
  es <- eigen(av)
  Ystdsample <- standardise_specifiedmultiplicity(Ysample, mult = c(2, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen(newav)
  expect_equal(newes$values[1], newes$values[2])
  expect_equal(newes$values[1], 2.5, tolerance = 0.3)
  expect_equal(newes$values[3], es$values[3])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(2,1)), 0)
  expect_error(expect_equal(stat_specifiedmultiplicity(Ysample, mult = c(2,1)), 0))
})

test_that("stat_specifiedmultiplicity() is zero for standarised sample, dim 7", {
  set.seed(13131)
  Ysample <- rsymm(50, 7)
  Ysample <- lapply(Ysample, `+`, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  av <- mmean(Ysample)
  es <- eigen(av)
  Ystdsample <- standardise_specifiedmultiplicity(Ysample, mult = c(3, 2, 1, 1))

  #check new average has correct properties
  newav <- mmean(Ystdsample)
  newes <- eigen(newav)
  expect_equal(newes$values[1:3], rep(3, 3), tolerance = 0.1)
  expect_equal(newes$values[4:5], rep(2, 2), tolerance = 0.1)
  expect_equal(newes$values[6], es$values[6])
  expect_equal(newes$values[7], es$values[7])
  expect_equal(newav %*% es$vectors %*% diag(1/newes$values), es$vectors)

  expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(3, 2, 1, 1)), 0)
  expect_error(expect_equal(stat_specifiedmultiplicity(Ystdsample, mult = c(3, 2, 2)), 0))
  expect_error(expect_equal(stat_specifiedmultiplicity(Ysample, mult = c(3, 2, 1, 1)), 0))
})

test_that("stat_specifiedevals() doesn't reject for simulation of single sample from null, and rejects otherwise", {
  set.seed(1331)
  Ysample <- rsymm(50, 7)
  Ysample <- lapply(Ysample, `+`, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)))
  res <- test_specifiedmultiplicity(Ysample, mult = c(3,2,1,1), 100)
  expect_gt(res$pval, 0.2)
})

})
