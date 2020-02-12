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
  estDL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                          lambda = 0.1, tree = tree,
                          alpha = 1, method = "desparsifiedlasso")
})

######################################################################
## Do not run this interactively. Insted, use Run Tests in RStudio. ##
######################################################################

test_that("outputs do not change over time", {
  # Shooting
  expect_known_value(estS$zscores_est, "previous_outputs/estS_zscore",
                     update = FALSE)
  expect_known_value(estS$shift_est, "previous_outputs/estS_shift",
                     update = FALSE)
  # L-BFGS-B
  expect_known_value(estL$zscores_est, "previous_outputs/estL_zscore",
                     update = FALSE)
  expect_known_value(estL$shift_est, "previous_outputs/estL_shift",
                     update = FALSE)
  # Scaled Lasso
  expect_known_value(estSL$zscores_est, "previous_outputs/estSL_zscore",
                     update = FALSE)
  expect_known_value(estSL$shift_est, "previous_outputs/estSL_shift",
                     update = FALSE)
  expect_known_value(estSL$optim_info$sigma_scaledlasso,
                     "previous_outputs/estSL_sigma", update = FALSE)
  # Desparsified Lasso
  expect_known_value(estDL$zscores_est, "previous_outputs/estDL_zscore",
                     update = FALSE)
  expect_known_value(estDL$shift_est, "previous_outputs/estDL_shift",
                     update = FALSE)
})


