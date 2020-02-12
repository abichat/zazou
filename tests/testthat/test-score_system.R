context("Score system")

m <- 10
n <- 100

X <- matrix(rnorm(m*n), nrow = n, ncol = m)

withr::with_preserve_seed({
  set.seed(42)
  scosys <- calculate_Z(X)
  set.seed(42)
  scosyshdi <- hdi:::calculate.Z(x = X, parallel = FALSE, ncores = 1,
                                 verbose = FALSE, Z = NULL, do.ZnZ = FALSE)
})

test_that("Score system has the right dimension", {
  expect_equal(dim(X), dim(scosys$Z))
})

test_that("Coherent with hdi code", {
  expect_equal(scosys$Z, scosyshdi$Z)
})

