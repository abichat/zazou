context("Matrices from trees")

library(ape)

tree <- read.tree(text = "((A:0.7, B:0.7):1.3, ((C:0.2, D:0.2):0.8, E:1):1);")
h <- 2

expected_incidence <- matrix(c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                     0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0, 0, 1), nrow = 5)


test_that("the incidence matrix is correct", {
  incidence <- incidence_matrix(tree)
  expect_equivalent(incidence + 0, expected_incidence)
  expect_equal(rownames(incidence), tree$tip.label)
})

test_that("the postorder is required for incidence matrix ", {
  expect_error(recursion_up(tree))
})


test_that("the covariance matrix is correct", {
  a <- rpois(1, lambda = 3) + 1
  covar <- covariance_matrix(tree, alpha = a)
  expect_equal(covar, t(covar))
  expect_equivalent(diag(covar), rep(1, length(tree$tip.label)))
  expect_equal(covar[1, 2],
               (exp(-a * 2 * 0.7) - exp(- 2 * a * h)) / (1 - exp(- 2 * a * h)),
               tolerance = 1e-5)
})
