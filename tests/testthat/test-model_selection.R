context("Hyperparameter selection")

data(chlamydiae)


pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)

# First tests

grid <- sample(c(1, 3))


estS <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = grid, tree = tree,
                        alpha = grid, method = "shooting")

estS2 <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         tree = tree, alpha = grid, method = "shooting")

estS3 <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         tree = tree, alpha = 1, lambda = 2,
                         method = "shooting", allow_positive = TRUE,
                         unknow = 3)

test_that("a selection is done or not", {
  expect_equal(ncol(estS$optim_info$bic_selection), 6)
  expect_equal(nrow(estS$optim_info$bic_selection), length(grid) ^ 2)
  expect_equal(nrow(estS2$optim_info$bic_selection),
               length(grid) * formals(lambda_grid)$n_lambda)
  expect_null(nrow(estS3$optim_info$bic_selection))
  expect_equal(estS$optim_info$criterion, "bic")
  expect_true(grepl("with model selection", estS$method))
  expect_true(grepl("with model selection", estS2$method))
  expect_false(grepl("with model selection", estS3$method))
})

test_that("supplementary arguments are tracked", {
  expect_length(estS$optim_info$supp_arg, 0)
  expect_length(estS3$optim_info$supp_arg, 2)
  expect_equal(names(estS3$optim_info$supp_arg),
               c("allow_positive", "unknow"))
})

# Second tests

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
  expect_equal(estL$shifts_est, estL_best$shifts_est)
  expect_true(estL_best$bic <= estL_notbest$bic)
  expect_lte(estL$bic, min(df_selection$bic))
})

est_pbic <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                            lambda = grid, tree = tree, criterion = "pbic",
                            alpha = grid, method = "shooting")

test_that("selection on pbic criterion is OK", {
  expect_true(grepl("with model selection", est_pbic$method))
  expect_equal(est_pbic$pbic, min(est_pbic$optim_info$bic_selection$pbic))
  expect_equal(est_pbic$optim_info$criterion, "pbic")
})

test_that("the shifts inside `bic_selection` are correct", {
  expect_equal(df_selection[df_selection$alpha == estL$alpha &
                              df_selection$lambda == estL$lambda,
                            "shifts_est"][[1]],
               estL$shifts_est)
  expect_equal(df_selection[df_selection$alpha == alpha_notbest &
                              df_selection$lambda == lambda_notbest, 6][[1]],
               estL_notbest$shifts_est)
})

all_est <- extract_models(estL)

r <- sample(seq_len(length(grid) ^ 2), size = 1)
est_r <- all_est[[r]]
est_r_fs <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                            lambda = est_r$lambda, tree = tree,
                            alpha = est_r$alpha, method = "L-BFGS-B")

test_that("extraction works correctly", {
  expect_length(all_est, length(grid) ^ 2)
  expect_true(grepl(", part of model selection", est_r$method))
  expect_equal(est_r$zscores_est, est_r_fs$zscores_est)
  expect_equal(est_r$objective_value, est_r_fs$objective_value)
  expect_equal(unname(est_r$optim_info$better_parameters["better_alpha"]),
               estL$alpha)
  expect_equal(unname(est_r$optim_info$better_parameters["better_lambda"]),
               estL$lambda)
  expect_warning(extract_models(estL_best))
})
