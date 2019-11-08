context("Unidirectionnal solver")


test_that("1-D solution is always negative", {
  expect_true(solve_univariate(rnorm(10), rnorm(10)) <= 0)
  expect_true(solve_univariate(rnorm(10), rnorm(10)) <= 0)
  expect_true(solve_univariate(rnorm(10), rnorm(10)) <= 0)
})


delta <- -10
n_change <- 5
n_keep <- 2

y <- c(rep(delta, n_change), rep(0, n_keep))
x <- c(rep(1, n_change), rep(0, n_keep))

test_that("1-D solution has the right boundary", {
  expect_equal(solve_univariate(y = y, x = x, lambda = 0), delta)
  expect_equal(solve_univariate(y = y, x = x, lambda = -crossprod(y, x)), 0)
  expect_true(solve_univariate(y, x, lambda = -crossprod(y, x) - 0.1) < 0)
})

test_that("1-D solution works when allowing positive values", {
  delta <- -10
  n_change <- 5
  n_keep <- 2

  y <- c(rep(delta, n_change), rep(0, n_keep))
  x <- c(rep(-1, n_change), rep(0, n_keep))

  expect_equal(solve_univariate(y = y, x = x, lambda = 0, TRUE), -delta)
  expect_equal(solve_univariate(y = y, x = x, lambda = crossprod(y, x), TRUE), 0)
  expect_gt(solve_univariate(y, x, lambda = crossprod(y, x) - 0.1, TRUE), 0)
})
