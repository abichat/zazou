context("Z-scores")

q95 <- qnorm(0.95)
q975 <- qnorm(0.975)
q025 <- qnorm(0.025)

test_that("unsigned z-scores are correct for remarkable values", {
  expect_equal(p2z(0.5), 0)
  expect_equal(p2z(0.05), q95)
  expect_equal(p2z(0.025), q975)
})

test_that("signed z-scores are correct for remarkable values", {
  expect_equal(p2z(0.05, e.sign = 1), q975)
  expect_equal(p2z(0.05, e.sign = -1), q025)
})

test_that("p2z is vectorized", {
  n <- 10
  expect_equal(p2z(rep(0.5, n)), rep(0, n))
})
