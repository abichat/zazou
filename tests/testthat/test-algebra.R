context("Algebra")

dim <- 26
M <- matrix(2, ncol = dim, nrow = dim) + diag(5 * dim, ncol = dim, nrow = dim)
R <- inverse_sqrt(M)

test_that("inverse_sqrt() is correct", {
  expect_equal(R[upper.tri(R)], rep(0, (dim * (dim-1))/2))
  expect_equal(t(R) %*% R %*% M, diag(1, ncol = dim, nrow = dim))
})


x0 <- rep(0, 5)
x1 <- 0:6
x2 <- c(-1, 5, 8)
x3 <- rnorm(100)
x3[sample(100, 10)] <- 0
x4 <- c(TRUE, FALSE)

test_that("norm0 is correct", {
  expect_equal(norm0(x0), 0)
  expect_equal(norm0(x1), 6)
  expect_equal(norm0(x1, reverse = TRUE), 1)
  expect_equal(norm0(x2), 3)
  expect_equal(norm0(x3), 90)
  expect_equal(norm0(x3, reverse = TRUE), 10)
  expect_error(norm0(x4))
})

