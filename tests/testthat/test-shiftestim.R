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
# estDL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
#                          lambda = 1, tree = tree,
#                          alpha = 1, method = "desparsifiedlasso")

estimations <- list(estS, estL, estSL)
estR <- estimations[[sample(x = length(estimations), size = 1)]]
pR <- plot(estR)

expect_scalnum <- function(x){
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
}

# General

test_that("Shiftestim class is correct", {
  expect_error(as_shiftestim(tree),
               "'listopt' must be the output of a optimisation function.")
  expect_equal(class(estR), "shiftestim")
  expect_is(estR$zscores_obs, "numeric")
  expect_is(estR$zscores_est, "numeric")
  expect_is(estR$shift_est, "data.frame")
  expect_is(estR$method, "character")
  expect_is(estR$optim_info, "list")
  expect_is(estR$optim_info$supp_arg, "list")
  expect_is(estR$tree, "phylo")
  expect_is(estR$is_bin, "logical")
  expect_scalnum(estR$lambda)
  expect_scalnum(estR$alpha)
  expect_scalnum(estR$sigma)
  expect_scalnum(estR$objective_value)
  expect_scalnum(estR$bic)
  expect_scalnum(estR$pbic)
  expect_scalnum(estR$pars_score)
  expect_is(pR, "ggtree")
  expect_is(pR, "gg")
  expect_is(pR, "ggplot")
})


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
})

