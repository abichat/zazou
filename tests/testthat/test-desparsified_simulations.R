context("Example desparsified with a tree")

n <- 15
tree <- ape::rtree(n)
tree <- force_ultrametric(tree)
nplusm <- length(tree$edge.length)
zscores <- simu_zscores(tree, 1, shifts = NULL, Nshifts = 3)

# High level

ESD <- estimate_shifts(Delta0 = rep(0, nplusm), zscores = zscores,
                       tree = tree, alpha = 1, lambda = 100,
                       method = "desparsified", alpha_conf = 0.01)


test_that("ESD has its specific components / dimensions", {
  expect_equal(ESD$method, "desparsified lasso")
  expect_equal(ESD$optim_info$alpha_confint, 0.01)
  expect_equal(ESD$optim_info$supp_arg$alpha_conf, 0.01)
  expect_true(is.na(ESD$objective_value))
  expect_equal(ncol(ESD$shift_est), 3)
  expect_is(ESD$optim_info$covariance_noise_matrix, "matrix")
})


mat_covar <- covariance_matrix(tree, alpha =  1)
mat_incidence <- incidence_matrix(tree)
R <- inverse_sqrt(mat_covar)
Y <- R %*% zscores
X <- R %*% mat_incidence
scosys <- calculate_Z(X = X)
tau <- noise_factor(X = X, score_system = scosys)
V <- covariance_noise_matrix(X = X, score_system = scosys)
j <- sample(nplusm, 1)
k <- sample(nplusm, 1)

scla <- solve_scaled_lasso(y = Y, X = X, beta0 = rep(0, nplusm),
                           constraint_type = "beta")
beta <- update_beta(X = X, y = Y, beta_init = scla$par$estimate,
                    score_system = scosys)
hci <- size_half_confint_shifts(noise_factor = tau,
                                hsigma = scla$sigma_scaledlasso)$half_size


test_that("intermediary outputs are correct", {
  ## Deterministic
  # score system
  expect_equal(dim(scosys), dim(X))
  # noise factor
  expect_length(tau, nplusm)
  # covariance noise matrix
  expect_equal(V, ESD$optim_info$covariance_noise_matrix)
  expect_equal(dim(V), c(nplusm, nplusm))
  expect_equal(V, t(V))
  expect_equal(V[j, k],
               sum(scosys[, j] * scosys[, k]) /
                 (abs(sum(scosys[, j] * X[, j]) * sum(scosys[, k] * X[, k]))))
  ## Deterministic
  # scaled lasso
  expect_length(scla, 6)
  expect_length(scla$par$estimate, nplusm)
  expect_length(scla$sigma_scaledlasso, 1)
  expect_named(scla, c("par", "sigma_scaledlasso", "method", "value",
                       "iterations", "last_progress"))
  # Updated beta
  expect_length(beta, nplusm)
  # Confidence interval for shifts
  expect_length(hci, nplusm)
})





