context("Column-wise inverse matrix")


test_that("estimation of polynom coefficients is correct for random 2*2", {
  m <- rnorm(2)
  a <- rnorm(4)
  A <- matrix(a, nrow = 2)
  coef_m1 <- coef_p2(A = A, m = m, ind = 1)
  coef_m2 <- coef_p2(A = A, m = m, ind = 2)
  expect_equal(coef_m1, list(a = a[1],
                             b = (a[2] + a[3]) * m[2],
                             c = a[4] * m[2]^2))
  expect_equal(coef_m2, list(a = a[4],
                             b = (a[2] + a[3]) * m[1],
                             c = a[1] * m[1]^2))
})

test_that("estimation of polynom coefficients is correct for 3*3", {
  m <- 1:3
  A <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
  coefs <- coef_p2(A = A, m = m, ind = 3)
  expect_equal(coefs, list(a = 9, b = 38, c = 33))
})

dim <- 10
gamma <- 2 * sqrt(log(dim)/dim)
# gamma <- 1
X <- matrix(rnorm(dim^2), ncol = dim) / dim
A <- t(X) %*% X
# M <- try(solve_colwiseinverse(A, gamma))

# while(!is.matrix(M)){
#   M <- try(solve_colwiseinverse(A, gamma))
# }
#
# test_that("global_constrains are repected", {
#   for(i in seq_len(dim)){
#     e <- rep(0, dim)
#     e[i] <- 1
#     expect_lte(max(A %*% M[, i] - e), gamma)
#   }
# })
