test_that("vec() works", {
  expect_equal(vec(matrix(1:12, byrow = FALSE, nrow = 3)),
               1:12)
})
test_that("invvec() works", {
  expect_equal(invvec(1:12, nrow = 3),
               matrix(1:12, byrow = FALSE, nrow = 3))
})
