test_that("dup() gives matrix that works on vec() and vech() results", {
  m <- matrix(c(1,2,3,2,5,6,3,6,9), byrow = FALSE, nrow = 3)
  expect_equal(drop(dup(3) %*% vech(m)), vec(m))
})

test_that("dup() memoisation works", {
  dup1 <- dup_direct(10)
  expect_equal(dup(10), dup1)
  time1 <- system.time(dup_direct(50))
  expect_lt(system.time(dup(50))["user.self"], time1["user.self"])
})
