context("Tree helpers")

library(ape)

tree1 <- rcoal(7)
tree2 <- rcoal(16)
tree3 <- rtree(42)
not_term_branches <- which(! tree3$edge[, 2] %in% seq_along(tree3$tip.label))

test_that("tree_height() is correct", {
  expect_equal(tree_height(tree1), max(cophenetic(tree1))/2)
  expect_equal(tree_height(tree2), max(cophenetic(tree2))/2)
  expect_error(tree_height(tree3))
})


test_that("fitch works", {
  tree <- read.tree(text = "(((t4:0.07418493367,(t6:0.0466889315,t5:0.0466889315):0.02749600216):0.6619201357,t2:0.7361050694):0.4046613234,(t7:0.3876200329,(t1:0.09560661687,t3:0.09560661687):0.2920134161):0.7531463598);")
  states <- c("t1" = "A",
              "t2" = "B",
              "t3" = "A",
              "t4" = "C",
              "t5" = "C",
              "t6" = "C",
              "t7" = "A"
              )
  expect_equal(fitch(tree, states), 2)
  states["t2"] <- "C"
  expect_equal(fitch(tree, states), 1)
  expect_warning(fitch(tree, unname(states)),
                 "State vector is unnamed, assuming same order as tip labels",
                 fixed = TRUE)
})

test_that("force_ultrametric() is correct", {
  expect_true(is.ultrametric(force_ultrametric(tree3)))
  expect_equal(force_ultrametric(tree3)$tip.label, tree3$tip.label)
  expect_equal(force_ultrametric(tree3)$edges, tree3$edges)
  expect_equal(force_ultrametric(tree3)$edge.length[not_term_branches],
               tree3$edge.length[not_term_branches])
})
