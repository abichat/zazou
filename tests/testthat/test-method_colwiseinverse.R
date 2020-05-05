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
M <- try(solve_colwiseinverse(A, gamma))

while(!is.matrix(M)){
  M <- try(solve_colwiseinverse(A, gamma))
}

test_that("global_constrains are repected", {
  for(i in seq_len(dim)){
    e <- rep(0, dim)
    e[i] <- 1
    expect_lte(max(A %*% M[, i] - e), gamma)
  }
})

test_that("Computation of m works for easy cases", {
  dim <- 10
  A <- diag(1, dim)
  for(i in seq_len(dim)){
    e <- rep(0, dim)
    e[i] <- 1
    ## Small gamma
    expect_identical(
      fast_solve_colwiseinverse_col(col = i, A = A, gamma = 0.1, m0 = e),
      matrix((1 - gamma) * e, ncol = 1)
    )
    ## Large gamma
    expect_identical(
      fast_solve_colwiseinverse_col(col = i, A = A, gamma = 1, m0 = e),
      matrix(0, nrow = dim, ncol = 1)
    )
  }
})

test_that("fast_solve_colwiseinverse_col() returns a warning/error
          when the initial vector is not feasible.", {
  dim <- 10
  A <- diag(1, dim)
  gamma <- 0.1
  expect_warning(try(fast_solve_colwiseinverse_col(col = 1, A = A, gamma = 0.1, m0 = rep(1, dim)), silent = TRUE),
                 "The starting point is not in the feasible set. Updates may be meaningless.",
                 fixed = TRUE)
  expect_error(suppressWarnings(fast_solve_colwiseinverse_col(col = 1, A = A, gamma = 0.1, m0 = rep(1, dim))),
               "No cell could be updated for the first time in column 1.",
               fixed = TRUE)
})

test_that("global_constrains are respected", {
  A <- toeplitz(x = 0.5^seq(0, 9))
  # B <- toeplitz(x = c(4/3, -2/3, rep(0, 8))) + diag(rep(c(0, 1/3, 0), times = c(1, 8, 1))) ## = solve(A)
  gamma <- 0.01
  m <- fast_solve_colwiseinverse_col(col = 1, A = A, gamma = gamma)
  e1 <- c(1, rep(0, 9))
  expect_lte(max(abs(A %*% m - e1)), gamma + sqrt(.Machine$double.eps))
})


test_that("giving A or svd(A) throws the same output", {
  A <- matrix(c(2, 4, 4, 9), ncol = 2)
  svdA <- svd(A)
  r1 <- withr::with_seed(2, fast_solve_colwiseinverse_col(col = 1, A = A,
                                                          gamma = 0.02))
  r2 <- withr::with_seed(2, fast_solve_colwiseinverse_col(col = 1, svdA = svdA,
                                                          gamma = 0.02))
  expect_equal(r1, r2)
})

test_that("find_feasible() returns a vector with correct dimensions", {
  B <- diag(1, 3)
  expect_length(object = find_feasible(B, 1, 0), n = ncol(B))
  B <- matrix(c(1, 0, 0, 1, 0, 0), ncol = 3)
  expect_length(object = find_feasible(B, 1, 0), n = ncol(B))
})

test_that("find_feasible() returns correct solution on easy cases", {
  B <- diag(1, 3)
  for (i in 1:ncol(B)) {
    expect_equal(object   = find_feasible(B, col = i, 0),
                 expected = rep(c(0, 1, 0), times = c(i - 1, 1, ncol(B) - i)))
  }
})

test_that("find_feasible() returns feasible solution", {
  B <- matrix(c(-1L, 0L, 1L, 1L, -1L, 3L, 1L, 1L, 0L), nrow = 3L) ## invertible matrix
  gamma <- 0.1
  for (i in 1:ncol(B)) {
    m <- find_feasible(B, i, gamma)
    e_i <- rep(c(0, 1, 0), times = c(i - 1, 1, ncol(B) - i))
    expect_lte(max(abs(B %*% m - e_i)), gamma + 1e-14)
  }
})

test_that("find_feasible() throws a warning when it does not find a feasible solution", {
  B <- matrix(c(1, 1, 0, 0), 2)
  gamma <- 0.4 ## Problem has no solution for gamma < 0.5
  expect_warning(find_feasible(B, 1, gamma = 0.4),
                 "No feasible solution found after 1e+05 iterations. Aborting",
                 fixed = TRUE)
  expect_warning(find_feasible(B, 1, gamma = 0.4, max_it = 10),
                 "No feasible solution found after 10 iterations. Aborting",
                 fixed = TRUE)
})





