context("Desparsified internal")

X <- matrix(1:12, ncol = 3)
Z <- matrix(rpois(12, 3), ncol = 3)
tau <- noise_factor_scoresystem(X, Z)

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






