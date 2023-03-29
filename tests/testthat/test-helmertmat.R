test_that("Helmert matrix has the correct properties", {
  # first row parallel to 11111
  # triangle of zeros
  # triangle of constant by row
  # orthonormal matrix
  n <- 5
  m <- helmert(n)
  expect_equal(m[1, ], rep(1/sqrt(n), n))
  expect_equal(m[-1, ][upper.tri(m, diag = FALSE)[-1, ]], rep(0, (n-2)*(n-1)/2))
  for (row in 2:n) {
    expect_equal(m[row, 1:(row-1)], rep(1, row-1)/sqrt(row-1 + (row-1)^2))
    expect_equal(m[row, row], -(row-1)/sqrt(row-1 + (row-1)^2))
  }

  expect_equal(m %*% t(m), diag(rep(1, n)))
  expect_equal(t(m) %*% m, diag(rep(1, n)))
})
