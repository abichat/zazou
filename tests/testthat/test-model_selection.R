test_that("lambda_grid works", {
  x <- matrix(c(1, 1, -1, 1), nrow = 2, ncol = 2)
  y <- c(-1/2, 1/2)
  lambda_max <- 1
  expect_equal(lambda_grid(x, y, n_lambda = 6, min_ratio = 1e-5),
               10^(-(0:5)))
  expect_equal(lambda_grid(x, y, n_lambda = 10, min_ratio = 1e-9),
               10^(-(0:9)))
})
