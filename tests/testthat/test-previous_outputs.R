context("Previous outputs")

data(chlamydiae)

pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)


withr::with_preserve_seed({
  set.seed(42)
  estS <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                          lambda = c(1, 2), tree = tree,
                          alpha = c(0.1, 2), method = "shooting")
  estL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                          lambda = 1, tree = tree,
                          alpha = 1, method = "L-BFGS-B")
  estSL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                          lambda = 0.1, tree = tree,
                          alpha = 1, method = "scaledlasso")
  estDL <- estimate_confint(estSL, alpha_conf = 0.05,
                            method = "desparsified")
})

######################################################################
## Do not run this interactively. Insted, use Run Tests in RStudio. ##
######################################################################

test_that("outputs do not change over time", {
  # Shooting
  expect_known_value(estS$zscores_est, "previous_outputs/estS_zscore",
                     update = FALSE)
  expect_known_value(estS$shifts_est, "previous_outputs/estS_shift",
                     update = FALSE)
  expect_known_value(estS$optim_info$bic_selection,
                     "previous_outputs/estS_bicselection",
                     update = FALSE)
  # L-BFGS-B
  expect_known_value(estL$zscores_est, "previous_outputs/estL_zscore",
                     update = FALSE)
  expect_known_value(estL$shifts_est, "previous_outputs/estL_shift",
                     update = FALSE)
  # Scaled Lasso
  expect_known_value(estSL$zscores_est,
                     "previous_outputs/estSL_zscore", update = FALSE)
  expect_known_value(estSL$shifts_est, "previous_outputs/estSL_shift",
                     update = FALSE)
  expect_known_value(estSL$optim_info$sigma_scaledlasso,
                     "previous_outputs/estSL_sigma", update = FALSE)
  # Desparsified Lasso
  expect_known_value(estDL$zscores_est,
                     "previous_outputs/estDL_zscore", update = FALSE)
  expect_known_value(estDL$shifts_est, "previous_outputs/estDL_shift",
                     update = FALSE)
  expect_known_value(estDL$optim_info$covariance_noise_matrix,
                     "previous_outputs/estDL_conoma",
                     update = FALSE)
  expect_known_value(estDL$shiftestim$zscores_est,
                     "previous_outputs/estSL_zscore",
                     update = FALSE)
})


