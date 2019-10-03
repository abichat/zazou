context("Inverse of the square root")

dim <- 26
M <- matrix(2, ncol = dim, nrow = dim) + diag(5 * dim, ncol = dim, nrow = dim)
R <- inverse_sqrt(M)

test_that("inverse_sqrt() is correct", {
  expect_equal(R[upper.tri(R)], rep(0, (dim * (dim-1))/2))
  expect_equal(t(R) %*% R %*% M, diag(1, ncol = dim, nrow = dim))
})
