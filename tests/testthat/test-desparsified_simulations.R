context("Example desparsified with a tree")

n <- 15
tree <- ape::rtree(n)
tree <- force_ultrametric(tree)
nplusm <- length(tree$edge.length)
zscores <- simu_zscores(tree, 1, shifts = NULL, Nshifts = 3)


ESD <- estimate_shifts(Delta0 = rep(0, nplusm),
                       zscores = zscores, tree = tree, alpha = 1, lambda = 100,
                       method = "desparsified", alpha_conf = 0.01)

plot(ESD)

mat_covar <- covariance_matrix(tree, alpha =  1)
mat_incidence <- incidence_matrix(tree)
R <- inverse_sqrt(mat_covar)
Y <- R %*% zscores
X <- R %*% mat_incidence
scosys <- calculate_Z(X = X)


test_that("ESD has its specific components / dimensions", {
  expect_equal(ESD$method, "desparsified lasso")
  expect_equal(ESD$optim_info$alpha_confint, 0.01)
  expect_equal(ESD$optim_info$supp_arg$alpha_conf, 0.01)
  expect_true(is.na(ESD$objective_value))
  expect_equal(ncol(ESD$shift_est), 3)
  expect_is(ESD$optim_info$covariance_noise_matrix, "matrix")
})

j <- sample(nplusm, 1)
k <- sample(nplusm, 1)
V <- ESD$optim_info$covariance_noise_matrix

test_that("the covariance noise matrix is correct", {
  expect_equal(dim(V), c(nplusm, nplusm))
  expect_equal(V, t(V))
  expect_equal(V[j, k],
               sum(scosys[, j] * scosys[, k]) /
                 (abs(sum(scosys[, j] * X[, j]) * sum(scosys[, k] * X[, k]))))

})





