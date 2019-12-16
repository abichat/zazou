context("Example desparsified for dimensions")

m <- 10
n <- 100

X <- matrix(rnorm(m*n), nrow = n, ncol = m)
y <- - 3 * x[, 1] - 5 * x[, 2] + 4 * x[, 3] + rnorm(n)

scla <- scaled_lasso(y = y, X = x)

test_that("scaled_lasso() has correct dimensions", {
  expect_length(scla, 2)
  expect_length(scla$beta_init, m)
  expect_length(scla$hsigma, 1)
  expect_named(scla, c("beta_init", "hsigma"))
  expect_true(all(scla$beta_init <= 0))
})

scosys <- score_system(X = X, y = y, beta_init = scla$beta_init,
                       hsigma = scla$hsigma)

test_that("score_system() has correct dimensions", {
  expect_equal(dim(scosys), dim(X))
})

tau <- noise_factor(X = X, score_system = scosys)

test_that("noise_factor() has correct dimensions", {
  expect_length(tau, m)
})

beta <- beta(X = X, y = y, beta_init = scla$beta_init, score_system = scosys)

test_that("noise_factor() has correct dimensions", {
  expect_length(beta, m)
})

hci <- size_half_confint(noise_factor = tau, hsigma = scla$hsigma)

test_that("size_half_confint() has correct dimensions", {
  expect_length(hci, m)
})

data.frame(lower = beta - hci, estimate = beta,
           upper = beta + hci, signif = beta > abs(hci))
