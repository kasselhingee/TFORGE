test_that("rsymm_t has correct mean", {
  set.seed(3654)
  sigma <- diag(c(3,2,1,1,1,1))
  Ys <- rsymm_t(1E4, delta = matrix(1, nrow = 3, ncol = 3), df = 10, sigma = sigma)
  expect_lt(max(abs(colMeans(Ys) - 1)), 2E-1)
  expect_equal(cov(Ys), sigma * 10/(10-2), 1E-1)
  # do.call(rbind, lapply(Ys, vech, name = TRUE)) |>
  #   data.frame() |>
  #   tidyr::pivot_longer(everything()) |>
  #   ggplot(aes(sample = value)) +
  #   facet_wrap(vars(name)) +
  #   geom_qq(pch = "+") +
  #   geom_qq_line(lty = "dashed", col = "blue")
})