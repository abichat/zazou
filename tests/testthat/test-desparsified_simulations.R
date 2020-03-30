context("Example desparsified with a tree")

n <- 15
tree <- ape::rtree(n)
tree <- force_ultrametric(tree)
nplusm <- length(tree$edge.length)
zscores <- simu_zscores(tree, 1, shifts = NULL, Nshifts = 3)

est_scaled <- estimate_shifts(Delta0 = rep(0, nplusm), zscores = zscores,
                              tree = tree, alpha = 1, lambda = 0.1,
                              method = "scaledlasso")


withr::with_preserve_seed({
  set.seed(42)
  ESD <- estimate_confint(est_scaled, alpha_conf = 0.05,
                          method = c("desparsified"))
})

withr::with_preserve_seed({
  set.seed(42)
  ESD002 <- estimate_confint(est_scaled, alpha_conf = 0.02,
                             method = c("desparsified"))
})


test_that("ESD has its specific components / dimensions", {
  expect_equal(ESD$method, "desparsified lasso")
  expect_equal(ESD$alpha_conf, 0.05)
  expect_equal(ncol(ESD$shifts_est), 3)
  expect_is(ESD$optim_info$covariance_noise_matrix, "matrix")
})


# mat_covar <- covariance_matrix(tree, alpha =  1)
# mat_incidence <- incidence_matrix(tree)
# R <- inverse_sqrt(mat_covar)
# Y <- R %*% zscores
# X <- R %*% mat_incidence
# scosys <- calculate_Z(X = X)
# tau <- noise_factor(X = X, score_system = scosys)
# V <- covariance_noise_matrix(X = X, score_system = scosys)
# j <- sample(nplusm, 1)
# k <- sample(nplusm, 1)
#
# scla <- solve_scaled_lasso(y = Y, X = X, beta0 = rep(0, nplusm),
#                            constraint_type = "beta")
# beta <- update_beta(X = X, y = Y, beta_init = scla$par$estimate,
#                     score_system = scosys)
# hci <- size_half_confint_shifts(noise_factor = tau,
#                                 hsigma = scla$sigma_scaledlasso)$half_size
# tau <- noise_factor(X, scosys)
#
#
# test_that("intermediary outputs are correct", {
#   ## Deterministic
#   # score system
#   expect_equal(dim(scosys), dim(X))
#   # noise factor
#   expect_length(tau, nplusm)
#   # covariance noise matrix
#   expect_equal(V, ESD$optim_info$covariance_noise_matrix)
#   expect_equal(dim(V), c(nplusm, nplusm))
#   expect_equal(V, t(V))
#   expect_equal(V[j, k],
#                sum(scosys[, j] * scosys[, k]) /
#                  (abs(sum(scosys[, j] * X[, j]) * sum(scosys[, k] * X[, k]))))
#   # Noise factor
#   expect_length(tau, nplusm)
#   expect_equal(tau, ESD$optim_info$noise_factor)
#   ## Non-deterministic
#   # scaled lasso
#   expect_length(scla, 6)
#   expect_length(scla$par$estimate, nplusm)
#   expect_length(scla$sigma_scaledlasso, 1)
#   expect_named(scla, c("par", "sigma_scaledlasso", "method", "value",
#                        "iterations", "last_progress"))
#   # Updated beta
#   expect_length(beta, nplusm)
#   # Confidence interval for shifts
#   expect_length(hci, nplusm)
# })


test_that("changing confindence interval works", {
  ESD5to2 <- update_confint(ESD, alpha_confint = 0.02)
  expect_equal(ESD5to2$alpha_conf, 0.02)
  expect_equal(ESD5to2$zscores_est$lower, ESD002$zscores_est$lower)
  expect_equal(ESD5to2$shifts_est$lower, ESD002$shifts_est$lower)
})


# ESD09 <- update_confint(ESD, alpha_confint = 0.9)
# ind <- which(ESD09$zscores_est$lower * ESD09$zscores_est$upper > 0)
#
# test_that("extraction works", {
#   expect_equal(extract_significant_leaves(ESD09, side = "both"),
#                ESD09$zscores_est$leaf[ind])
# })



