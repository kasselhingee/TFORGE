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



