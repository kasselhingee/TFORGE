test_that("vech() works", {
  expect_equal(vech(matrix(1:9, byrow = FALSE, nrow = 3)),
               c(1:3, 5:6, 9))
  expect_equal(
    names(vech(matrix(1:9, byrow = FALSE, nrow = 3), name = TRUE)),
    c("e11", "e21", "e31", "e22", "e32", "e33"))
})

test_that("inv_vech() works", {
  expect_equal(inv_vech(c(1:3, 5:6, 9)),
               matrix(c(1,2,3,2,5,6,3,6,9), byrow = FALSE, nrow = 3))
})
