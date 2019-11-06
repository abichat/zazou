context("Hyperparameter selection")

data(chlamydiae)


pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)

grid <- c(0.1, 1, 5, 7)

estS <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = grid, tree = tree,
                        alpha = grid, method = "shooting")

test_that("a selection is done", {
  expect_equal(nrow(estS$optim_info$bic_selection), length(grid)^2)
  expect_true(grepl("with model selection", estS$method))
})


estL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = grid, tree = tree,
                        alpha = grid, method = "L-BFGS-B")

estL_best <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                             lambda = estL$lambda, tree = tree,
                             alpha = estL$alpha, method = "L-BFGS-B")

estL_notbest <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                                lambda = sample(setdiff(grid, estL$lambda), 1),
                                alpha = sample(setdiff(grid, estL$alpha), 1),
                                tree = tree, method = "L-BFGS-B")

test_that("the choosen model is the best one", {
  expect_true(estL_best$bic <= estL_notbest$bic)
  expect_equal(estL$shift_est, estL_best$shift_est)
})




