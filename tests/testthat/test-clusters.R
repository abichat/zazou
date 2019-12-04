context("Clusters")

library(ape)
tree <- rtree(10)
N <- 4
cl_mono <- create_clusters(tree, N)
cl_para <- create_clusters(tree, N, "paraphyletic")
cl_unif <- create_clusters(tree, N, "uniform")

test_that("create_clusters() has the correct format", {
  expect_is(cl_mono, "integer")
  expect_is(cl_para, "integer")
  expect_is(cl_unif, "integer")
  expect_length(cl_mono, length(tree$tip.label))
  expect_length(cl_para, length(tree$tip.label))
  expect_length(cl_unif, length(tree$tip.label))
  expect_equal(sort(names(cl_mono)), sort(tree$tip.label))
  expect_equal(sort(names(cl_para)), sort(tree$tip.label))
  expect_equal(sort(names(cl_unif)), sort(tree$tip.label))
  expect_equal(length(unique(cl_mono)), N)
  expect_equal(length(unique(cl_para)), N)
  expect_equal(length(unique(cl_unif)), N)
  expect_equivalent(create_clusters(tree, 10), 1:10)
  expect_equivalent(create_clusters(tree, 1), rep(1, 10))
  expect_equal(length(unique(create_clusters(tree, 9, "uniform"))), 9)
})

test_that("create_clusters() throws errors", {
  expect_error(create_clusters(tree, 11))
  expect_error(create_clusters(tree, 0))
  expect_error(create_clusters(tree, -1))
})


