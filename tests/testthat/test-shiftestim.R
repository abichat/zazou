context("Shiftestim class")


data(chlamydiae)


pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)
N_branch <- length(tree$edge.length)

# First tests

estS <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = c(1, 2), tree = tree,
                        alpha = 1, method = "shooting")
estL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                        lambda = 1, tree = tree,
                        alpha = 1, method = "L-BFGS-B")
estSL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         lambda = 1, tree = tree,
                         alpha = 1, method = "scaledlasso")
estDL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                         lambda = 1, tree = tree, alpha_conf = 0.01,
                         alpha = 1, method = "desparsifiedlasso")

expect_scalnum <- function(x){
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
}

expect_shiftestim <- function(x){
  px <- plot(x)
  expect_equal(class(x), "shiftestim")
  expect_is(x$zscores_obs, "numeric")
  expect_is(x$zscores_est, "numeric")
  expect_is(x$shift_est, "data.frame")
  expect_is(x$method, "character")
  expect_is(x$optim_info, "list")
  expect_is(x$optim_info$supp_arg, "list")
  expect_is(x$tree, "phylo")
  expect_is(x$is_bin, "logical")
  expect_scalnum(x$lambda)
  expect_scalnum(x$alpha)
  expect_scalnum(x$sigma)
  # expect_scalnum(x$objective_value) # Not available for DL for the moment
  expect_scalnum(x$bic)
  expect_scalnum(x$pbic)
  expect_scalnum(x$pars_score)
  expect_is(px, "ggtree")
  expect_is(px, "gg")
  expect_is(px, "ggplot")
}


# General

expect_error(as_shiftestim(tree),
             "'listopt' must be the output of a optimisation function.")

test_that("Shiftestim class is correct", {
  expect_shiftestim(estS)
  expect_shiftestim(estL)
  expect_shiftestim(estSL)
  expect_shiftestim(estDL)
})

# Methods specificities

test_that("shooting output is correct", {
  expect_equal(estS$method, "shooting with model selection")
  expect_equal(ncol(estS$shift_est), 1)
  expect_scalnum(estS$optim_info$last_progress)
  expect_scalnum(estS$optim_info$iterations)
})

test_that("L-BFGS-B output is correct", {
  expect_equal(estL$method, "L-BFGS-B")
  expect_equal(ncol(estL$shift_est), 1)
  expect_is(estL$optim_info$message, "character")
  expect_error(
    estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                    lambda = 1, tree = tree, constraint_type = "yhat",
                    alpha = 1, method = "L-BFGS-B"),
    "The constraint 'yhat' is not available for L-BFGS-B solving.")
})

test_that("shooting output is correct", {
  expect_equal(estSL$method, "scaled lasso")
  expect_equal(ncol(estSL$shift_est), 1)
  expect_scalnum(estSL$optim_info$last_progress)
  expect_scalnum(estSL$optim_info$iterations)
  expect_error(
    estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
                    tree = tree, alpha = 1, method = "scaledlasso"),
    "No model selection can be done with scaledlasso.")
})

test_that("desparsified output is correct", {
  expect_equal(estDL$method, "desparsified lasso")
  expect_equal(ncol(estDL$shift_est), 3)
  expect_equal(estDL$optim_info$alpha_confint, 0.01)
  expect_true(is.na(estDL$objective_value))
})

