context("Incidence matrix")

library(ape)

tree <- read.tree(text = "((A, B), ((C, D), E));")

expected <- matrix(c(1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                     0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0, 0, 1), nrow = 5)

test_that("the incidence matrix is correct", {
  expect_equal(unname(incidence_matrix(tree)) + 0, expected)
  expect_equal(rownames(incidence_matrix(tree)), tree$tip.label)
})

test_that("the postorder is required", {
  expect_error(recursion_up(tree))
})

