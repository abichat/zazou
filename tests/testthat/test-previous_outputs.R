context("Previous outputs")

data(chlamydiae)

pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)

withr::with_seed(42, {
  estS <- estimate_shifts(zscores = zsco_obs,
                          lambda = c(1, 2), tree = tree,
                          alpha = c(0.1, 2), method = "shooting")
})

withr::with_seed(42, {
  estL <- estimate_shifts(zscores = zsco_obs,
                          lambda = 1, tree = tree,
                          alpha = 1, method = "L-BFGS-B")
})

withr::with_seed(42, {
  estSL <- estimate_shifts(zscores = zsco_obs,
                          lambda = 0.1, tree = tree,
                          alpha = 1, method = "scaledlasso")
})

withr::with_seed(42, {
  estSS <- estimate_confint(estSL, alpha_conf = 0.05,
                            method = "scoresystem")
})

withr::with_seed(42, {
  estCI <- suppressWarnings(estimate_confint(estSL, alpha_conf = 0.05,
                                             method = "colwiseinverse"))
})

######################################################################
## Do not run this interactively. Insted, use Run Tests in RStudio. ##
######################################################################

test_that("Sooting outputs do not change over time", {
  expect_known_value(estS$zscores_est, "previous_outputs/estS_zscore",
                     update = FALSE)
  expect_known_value(estS$shifts_est, "previous_outputs/estS_shift",
                     update = FALSE)
  expect_known_value(estS$optim_info$bic_selection,
                     "previous_outputs/estS_bicselection",
                     update = FALSE)
})


test_that("L-BFGS-B outputs do not change over time", {
  expect_known_value(estL$zscores_est, "previous_outputs/estL_zscore",
                     update = FALSE)
  expect_known_value(estL$shifts_est, "previous_outputs/estL_shift",
                     update = FALSE)
})


test_that("Scaled Lasso outputs do not change over time", {
  expect_known_value(estSL$zscores_est,
                     "previous_outputs/estSL_zscore", update = FALSE)
  expect_known_value(estSL$shifts_est, "previous_outputs/estSL_shift",
                     update = FALSE)
  expect_known_value(estSL$optim_info$sigma_scaledlasso,
                     "previous_outputs/estSL_sigma", update = FALSE)
})


test_that("Score system outputs do not change over time", {
  expect_known_value(estSS$zscores_est,
                     "previous_outputs/estSS_zscore", update = FALSE)
  expect_known_value(estSS$shifts_est, "previous_outputs/estSS_shift",
                     update = FALSE)
  expect_known_value(estSS$covariance_noise_matrix,
                     "previous_outputs/estSS_conoma",
                     update = FALSE)
  expect_known_value(estSS$optim_info$scoresystem,
                     "previous_outputs/estSS_scoresystem",
                     update = FALSE)
  expect_known_value(estSS$shiftestim$zscores_est,
                     "previous_outputs/estSL_zscore",
                     update = FALSE)
})


test_that("Columnwise inverse outputs do not change over time", {
  expect_known_value(estCI$zscores_est,
                     "previous_outputs/estCI_zscore", update = FALSE)
  expect_known_value(estCI$shifts_est, "previous_outputs/estCI_shift",
                     update = FALSE)
  expect_known_value(estCI$covariance_noise_matrix,
                     "previous_outputs/estCI_conoma",
                     update = FALSE)
  expect_known_value(estCI$optim_info$colwiseinverse,
                     "previous_outputs/estCI_scoresystem",
                     update = FALSE)
  expect_known_value(estCI$shiftestim$zscores_est,
                     "previous_outputs/estSL_zscore",
                     update = FALSE)
})


