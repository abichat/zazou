context("Hierarchical correction")

data(chlamydiae)

pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
tree <- force_ultrametric(chlamydiae$tree)

withr::with_seed(42, {
  zsco_obs <- p2z(pval_obs)
  est1 <- estimate_shifts(zscores = zsco_obs,
                           lambda = c(0.1, 2), tree = tree,
                           alpha = c(0.1, 2), method = "scaled lasso")
  est2 <- estimate_confint(est1, alpha_conf = 0.05,
                           method = "scoresystem")
})

withr::with_seed(42, {
  pv <- smooth_pvalues(pvalues = pval_obs, tree = tree,
                        arg_shiftestim = list(lambda = c(0.1, 2),
                                              alpha = c(0.1, 2),
                                              method = "scaled lasso"))
})


test_that("wrapper works", {
  expect_equal(unname(pv), est2$zscores_est$pvalue)
  expect_equal(names(pv), est2$zscores_est$leaf)
})
