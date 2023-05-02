test_that("rsymm_t has correct mean", {
  set.seed(3654)
  Ys <- rsymm_t(100, delta = matrix(1, nrow = 3, ncol = 3), df = 10, sigma = diag(c(3,2,1,1,1,1)))
  arr <- abind::abind(Ys, along = 3)
  mn <- apply(arr, MARGIN = c(1, 2), mean) #close to 1 - yes :)
  expect_lt(max(abs(mn - 1)), 2E-1)
  # do.call(rbind, lapply(Ys, vech, name = TRUE)) |>
  #   data.frame() |>
  #   tidyr::pivot_longer(everything()) |>
  #   ggplot(aes(sample = value)) +
  #   facet_wrap(vars(name)) +
  #   geom_qq(pch = "+") +
  #   geom_qq_line(lty = "dashed", col = "blue")
})