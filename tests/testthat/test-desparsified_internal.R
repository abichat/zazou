context("Desparsified internal")

X <- matrix(1:12, ncol = 3)
Z <- matrix(rpois(12, 3), ncol = 3)
tau <- noise_factor(X, Z)

test_that("noise_factor() is correct", {
  expect_length(tau, 3)
  # expect_equal(noise_factor(X, X), rep(1, ncol(X)))
  # expect_equal(noise_factor(Z, Z), rep(1, ncol(Z)))
  expect_equal(tau[1],
               norm(Z[, 1, drop = FALSE], type = "2") / sum(Z[, 1] * X[, 1]))
  expect_equal(tau[2],
               norm(Z[, 2, drop = FALSE], type = "2") / sum(Z[, 2] * X[, 2]))
  expect_equal(tau[3],
               norm(Z[, 3, drop = FALSE], type = "2") / sum(Z[, 3] * X[, 3]))
})


shc <- size_half_confint(1:4, 1)

test_that("size_half_confint() is correct", {
  expect_equal(shc[4], shc[1] * 4)
  expect_equal(size_half_confint(runif(5), 1, 0), rep(Inf, 5))
  expect_equal(size_half_confint(runif(5), 1, 1), rep(0, 5))
  expect_error(size_half_confint(runif(5), 1:4, 0))
  expect_error(size_half_confint(runif(5), 1, 0:1))
})




