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
})

######################################################################
## Do not run this interactively. Insted, use Run Tests in RStudio. ##
######################################################################

test_that("outputs do not change over time", {
  expect_known_value(estS$zscores_est, "previous_outputs/estS_zscore",
                     update = FALSE)
  expect_known_value(estL$zscores_est, "previous_outputs/estL_zscore",
                     update = FALSE)
  expect_known_value(estS$shift_est, "previous_outputs/estS_shift",
                     update = FALSE)
  expect_known_value(estL$shift_est, "previous_outputs/estL_shift",
                     update = FALSE)
})


