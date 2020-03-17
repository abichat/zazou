context("Desparsified internal")

X <- matrix(1:12, ncol = 3)
Z <- matrix(rpois(12, 3), ncol = 3)
tau <- noise_factor(X, Z)

test_that("noise_factor() is correct", {
  expect_length(tau, 3)
  # expect_equal(tau[1],
  #              norm(Z[, 1, drop = FALSE], type = "2") / sum(Z[, 1] * X[, 1]))
  # expect_equal(tau[2],
  #              norm(Z[, 2, drop = FALSE], type = "2") / sum(Z[, 2] * X[, 2]))
  # expect_equal(tau[3],
  #              norm(Z[, 3, drop = FALSE], type = "2") / sum(Z[, 3] * X[, 3]))
  expect_equal(tau[1], norm(Z[, 1, drop = FALSE], type = "2") / nrow(X))
  expect_equal(tau[2], norm(Z[, 2, drop = FALSE], type = "2") / nrow(X))
  expect_equal(tau[3], norm(Z[, 3, drop = FALSE], type = "2") / nrow(X))
})


shc <- size_half_confint_shifts(1:4, 1)

test_that("size_half_confint() is correct", {
  expect_equal(shc$half_size[4], shc$half_size[1] * 4)
  expect_equal(size_half_confint_shifts(runif(5), 1, 0)$half_size, rep(Inf, 5))
  expect_equal(size_half_confint_shifts(runif(5), 1, 1)$half_size, rep(0, 5))
  expect_equal(size_half_confint_shifts(runif(5), 1, alpha = 0.65)$alpha, 0.65)
  expect_error(size_half_confint_shifts(runif(5), 1:4, 0))
  expect_error(size_half_confint_shifts(runif(5), 1, 0:1))
})





