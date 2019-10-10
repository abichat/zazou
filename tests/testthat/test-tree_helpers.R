context("Tree helpers")

library(ape)

tree1 <- rcoal(7)
tree2 <- rcoal(16)
tree3 <- rtree(4)

test_that("tree_height() is correct", {
  expect_equal(tree_height(tree1), max(cophenetic(tree1))/2)
  expect_equal(tree_height(tree2), max(cophenetic(tree2))/2)
  expect_error(tree_height(tree3))
})

