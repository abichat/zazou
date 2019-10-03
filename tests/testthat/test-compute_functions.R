context("Compute objective function")

n <- 3
m <- 4

Y <- 10*1:m
X <- matrix(1:(m*(n+m)), nrow = m, ncol = n+m)
lambda <- 5
Delta <- -(1:(n+m))

t(Y - X %*% Delta) %*% (Y - X %*% Delta) - lambda * sum(Delta)


test_that("objective function is correct", {
  obj <- compute_objective_function(Y, X, lambda)
  value <- obj(Delta)
  expect_is(obj, "function")
  expect_length(value, 1)
  expect_equal(value,
               t(Y - X %*% Delta) %*% (Y - X %*% Delta) - lambda * sum(Delta))
})
