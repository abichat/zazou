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

test_that("1-D solution works when lifting the constraints", {
  delta <- -10
  n_change <- 5
  n_keep <- 2

  y <- c(rep(delta, n_change), rep(0, n_keep))
  x <- c(rep(-1, n_change), rep(0, n_keep))

  expect_equal(solve_univariate(y = y, x = x, lambda = 0,
                                constraint_type = "none"), -delta)
  expect_equal(solve_univariate(y = y, x = x, lambda = crossprod(y, x),
                                constraint_type = "none"), 0)
  expect_gt(solve_univariate(y, x, lambda = crossprod(y, x) - 0.1,
                             constraint_type = "none"), 0)
})

test_that("1-D solution works when allowing positive values and forgoing z", {
  delta <- -10
  n_change <- 5
  n_keep <- 2

  y <- c(rep(delta, n_change), rep(0, n_keep))
  x <- c(rep(-1, n_change), rep(0, n_keep))

  expect_equal(solve_univariate(y = y, x = x, u = -1, v = -1,
                                lambda = 0, constraint_type = "yhat"),
               -delta)
  expect_equal(solve_univariate(y = y, x = x, u = -1, v = 2,
                                lambda = 0, constraint_type = "yhat"),
               1/2) # projection on [-Inf, 1/2]
  expect_equal(solve_univariate(y, x, u = -1, v = -1,
                                lambda = crossprod(y, x),
                                constraint_type = "yhat"), 0)
  expect_gt(solve_univariate(y, x, u = -1, v = -1,
                             lambda = crossprod(y, x) - 0.01,
                             constraint_type = "yhat"), 0)
})

# test_that("1-D solution fails when the constraint is not satisfiable", {
#   # Can't find an unfeasible constraint
#   expect_error(solve_univariate(y = y, x = x, u = 10, v = 0,
#                                 lambda = 0, constraint_type = "yhat"),
#           "The constraint is not feasible. Consider changing the constraint.")
# })

# test_that("1-D solution works when the constraint is satisfiable", {
#   y <- c(-1, -1, -1, -1)
#   z <- c(-2, -1, -2, -1)
#   x <- c(1, 0, 1, 0)
#   # z + x = y
#   # Redundant, no?
#   expect_equal(solve_univariate(y = y, x = x, z = z,
#                                 constraint_type = "beta"), 0)
#   expect_equal(solve_univariate(y = y, x = x, z = z,
#                                 constraint_type = "yhat"), 1)
#   expect_equal(solve_univariate(y = y, x = x, z = z,
#                                 constraint_type = "none"), 1)
# })
