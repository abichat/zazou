context("Shiftestim class")


data(chlamydiae)


pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
zsco_obs <- p2z(pval_obs)

tree <- force_ultrametric(chlamydiae$tree)

# First tests

estS <- estimate_shifts(zscores = zsco_obs,
                        lambda = c(1, 2), tree = tree,
                        alphaOU = 1, method = "lasso")
estL <- estimate_shifts(zscores = zsco_obs,
                        lambda = 1, tree = tree,
                        alphaOU = 1, method = "L-BFGS-B")
estSL <- estimate_shifts(zscores = zsco_obs,
                         lambda = 1, tree = tree,
                         alphaOU = 1, method = "scaled lasso")

expect_scalnum <- function(x){
  expect_is(x, "numeric")
  expect_equal(length(x), 1)
}

expect_shiftestim <- function(x){
  px <- plot(x)
  expect_equal(class(x), "shiftestim")
  expect_is(x$zscores_obs, "numeric")
  expect_is(x$zscores_est, "numeric")
  expect_is(x$shifts_est, "numeric")
  expect_is(x$method, "character")
  expect_is(x$optim_info, "list")
  expect_is(x$optim_info$supp_arg, "list")
  expect_is(x$tree, "phylo")
  expect_is(x$is_bin, "logical")
  expect_scalnum(x$lambda)
  expect_scalnum(x$alpha)
  expect_scalnum(x$sigmaOU)
  expect_scalnum(x$objective_value)
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
  # expect_shiftestim(estDL)
})

# Methods specificities

test_that("lasso output is correct", {
  expect_equal(estS$method, "lasso with model selection")
  expect_scalnum(estS$optim_info$last_progress)
  expect_scalnum(estS$optim_info$iterations)
})

test_that("L-BFGS-B output is correct", {
  expect_equal(estL$method, "L-BFGS-B")
  expect_is(estL$optim_info$message, "character")
  expect_error(
    estimate_shifts(zscores = zsco_obs,
                    lambda = 1, tree = tree, constraint_type = "yhat",
                    alpha = 1, method = "L-BFGS-B"),
    "The constraint 'yhat' is not available for L-BFGS-B solving.")
})

test_that("scaled lasso output is correct", {
  expect_equal(estSL$method, "scaled lasso")
  expect_scalnum(estSL$optim_info$last_progress)
  expect_scalnum(estSL$optim_info$iterations)
  expect_scalnum(estSL$optim_info$sigma_scaledlasso)
})

# Methods

test_that("pull_pvalues() works", {
  expect_equal(pull_pvalues(estS), pnorm(estS$zscores_est))
})

# test_that("scoresystem output is correct", {
#   # expect_equal(estDL$method, "scoresystem")
#   # expect_equal(ncol(estDL$shifts_est), 3)
#   # expect_equal(estDL$optim_info$alpha_confint, 0.01)
#   # expect_true(is.na(estDL$objective_value))
#   # Check confidence interval
#   # expect_true(all(estDL$shifts_est$lower < estDL$shifts_est$estimate))
#   # expect_true(all(estDL$shifts_est$upper > estDL$shifts_est$estimate))
#   # expect_equal(estDL$shifts_est$upper - estDL$shifts_est$estimate,
#   #              estDL$shifts_est$estimate - estDL$shifts_est$lower)
#   # expect_true(all(estDL$zscores_est$lower < estDL$zscores_est$estimate))
#   # expect_true(all(estDL$zscores_est$upper > estDL$zscores_est$estimate))
#   # expect_equal(estDL$zscores_est$upper - estDL$zscores_est$estimate,
#   #              estDL$zscores_est$estimate - estDL$zscores_est$lower)
# })

# test_that("warnings and errors for shiftestims functions are correct", {
#   expect_error(update_confint(tree), "x must be a 'shiftestim' object.")
#   expect_warning(update_confint(estL))
#   # expect_warning(update_confint(estL),
#   #              "There is no confindence interval for this method (L-BFGS-B).")
# })

