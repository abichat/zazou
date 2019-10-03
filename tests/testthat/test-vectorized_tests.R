context("Vectorized tests")

X <- alcohol$X
Y <- alcohol$Y

i <- 1
j <- 10
k <- 100
r <- sample(nrow(X), size = 1)


wilcoxon <- test_wilcoxon(X, Y)
kruskalwallis <- test_kruskalwallis(X, Y)
fisher <- test_fisher(X, Y)

test_that("formats are correct", {
  expect_is(wilcoxon, "list")
  expect_equal(length(wilcoxon), 2)
  expect_equal(names(wilcoxon$p.value), rownames(X))
  expect_equal(names(wilcoxon$e.sign), rownames(X))


  expect_is(kruskalwallis, "list")
  expect_equal(length(kruskalwallis), 1)
  expect_equal(names(kruskalwallis$p.value), rownames(X))

  expect_is(fisher, "list")
  expect_equal(length(fisher), 1)
  expect_equal(names(fisher$p.value), rownames(X))
})

test_that("wilcoxon is correct", {
  expect_equal(unname(wilcoxon$p.value[c(i, j, k, r)]),
               c(suppressWarnings(wilcox.test(X[i, ] ~ Y)$p.value),
                 suppressWarnings(wilcox.test(X[j, ] ~ Y)$p.value),
                 suppressWarnings(wilcox.test(X[k, ] ~ Y)$p.value),
                 suppressWarnings(wilcox.test(X[r, ] ~ Y)$p.value)))
})

test_that("kruskalwallis is correct", {
  expect_equal(unname(kruskalwallis$p.value[c(i, j, k, r)]),
               c(suppressWarnings(kruskal.test(X[i, ] ~ Y)$p.value),
                 suppressWarnings(kruskal.test(X[j, ] ~ Y)$p.value),
                 suppressWarnings(kruskal.test(X[k, ] ~ Y)$p.value),
                 suppressWarnings(kruskal.test(X[r, ] ~ Y)$p.value)))
})

test_that("fisher is correct", {
  expect_equal(unname(fisher$p.value[c(i, j, k, r)]),
               c(suppressWarnings(anova(lm(X[i, ] ~ Y))$`Pr(>F)`[1]),
                 suppressWarnings(anova(lm(X[j, ] ~ Y))$`Pr(>F)`[1]),
                 suppressWarnings(anova(lm(X[k, ] ~ Y))$`Pr(>F)`[1]),
                 suppressWarnings(anova(lm(X[r, ] ~ Y))$`Pr(>F)`[1])))
})


Y2 <- as.character(Y)
Y2[1:3] <- "Medium"

test_that("wilcoxon throws error with more than two conditions", {
  expect_error(test_wilcoxon(X, Y2))
})

