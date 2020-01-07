context("Example desparsified for dimensions")

m <- 10
n <- 100

set.seed(42)

X <- matrix(rnorm(m*n), nrow = n, ncol = m)
X <- cbind(1, X)
y <- -30 - 3 * X[, 1] - 5 * X[, 2] + 4 * X[, 3] + rnorm(n)
# X <- scale(X, center = TRUE, scale = FALSE)
# y <- y - mean(y)

Beta0 <- rep(0, ncol(X))
scla <- scaled_lasso(y = y, X = X, projected = FALSE)
scaled_lasso2(y = y, X = X, beta0 = rep(0, ncol(X)), lambda = 10^-1, use_constraint = FALSE, allow_positive = TRUE)
beta <- scaled_lasso2(y = y, X = X, beta0 = rep(0, ncol(X)), lambda = 10^-1, use_constraint = FALSE, allow_positive = TRUE)$beta
scaled_lasso2(y = y, X = X, beta0 = beta, lambda = 10^-1, use_constraint = TRUE, allow_positive = TRUE)
max(X %*% beta)
scalreg(X, y, lam0 = 10^-1)

test_that("scaled_lasso() has correct dimensions", {
  expect_length(scla, 2)
  expect_length(scla$beta_init, m)
  expect_length(scla$hsigma, 1)
  expect_named(scla, c("beta_init", "hsigma"))
  # expect_true(all(scla$beta_init <= 0))
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
           upper = beta + hci, signif = abs(beta) > hci)

bhat <- lasso.proj(x = X, y = y,
                   # betainit = scla$beta_init, sigma = scla$hsigma,
                   betainit = "scaled lasso",
                   return.Z = TRUE, standardize = FALSE)$bhat

(bhat - beta)
(bhat - beta) / beta
# plot(bhat, beta)
