context("Hyperparameter selection")

data(chlamydiae)


pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)

grid <- sample(c(1, 3))

estS <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = grid, tree = tree,
                        alpha = grid, method = "shooting")

estS2 <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         tree = tree, alpha = grid, method = "shooting")

estS3 <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         tree = tree, alpha = 1, lambda = 2,
                         method = "shooting")

test_that("a selection is done or not", {
  expect_equal(ncol(estS$optim_info$bic_selection), 4)
  expect_equal(nrow(estS$optim_info$bic_selection), length(grid) ^ 2)
  expect_equal(nrow(estS2$optim_info$bic_selection), length(grid) * 6)
  expect_null(nrow(estS3$optim_info$bic_selection))
  expect_true(grepl("with model selection", estS$method))
  expect_true(grepl("with model selection", estS2$method))
  expect_false(grepl("with model selection", estS3$method))
})

grid <- sample(c(0.9, 1, 1.1))

estL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = grid, tree = tree,
                        alpha = sample(grid), method = "L-BFGS-B")

df_selection <- estL$optim_info$bic_selection

estL_best <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                             lambda = estL$lambda, tree = tree,
                             alpha = estL$alpha, method = "L-BFGS-B")

alpha_notbest <- sample(setdiff(grid, estL$alpha), 1)
lambda_notbest <- sample(setdiff(grid, estL$lambda), 1)

estL_notbest <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                                alpha = alpha_notbest,
                                lambda = lambda_notbest,
                                tree = tree, method = "L-BFGS-B")

test_that("the choosen model is the best one", {
  expect_equal(estL$shift_est, estL_best$shift_est)
  expect_true(estL_best$bic <= estL_notbest$bic)
  expect_lte(estL$bic, min(df_selection$bic))
})

test_that("the shifts inside `bic_selection` are correct", {
  expect_equal(df_selection[df_selection$alpha == estL$alpha &
                              df_selection$lambda == estL$lambda,
                            "shift_est"][[1]],
               estL$shift_est)
  expect_equal(df_selection[df_selection$alpha == alpha_notbest &
                              df_selection$lambda == lambda_notbest, 4][[1]],
               estL_notbest$shift_est)
})



