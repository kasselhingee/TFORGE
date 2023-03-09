test_that("vech() works", {
  expect_equal(vech(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1:3, 5:6, 9))
})

test_that("invvech() works", {
  expect_equal(invvech(c(1:3, 5:6, 9)),
               matrix(c(1,2,3,2,5,6,3,6,9), byrow = FALSE, nrow = 3))
})
