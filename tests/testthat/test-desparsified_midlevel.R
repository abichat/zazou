context("Example desparsified for dimensions")

m <- 10
n <- 100

set.seed(42)

X <- matrix(rnorm(m*n), nrow = n, ncol = m)
# X <- cbind(1, X)
y <- - 3 * X[, 1] - 5 * X[, 2] + 4 * X[, 3] + rnorm(n)
# X <- scale(X, center = TRUE, scale = FALSE)
# y <- y - mean(y)

# Beta0 <- rep(0, ncol(X))
scalreg::scalreg(X, y)
(scla <- solve_scaled_lasso(y = y, X = X, beta0 = rep(0, ncol(X)),
                            constraint_type = "none"))
(scla <- solve_scaled_lasso(y = y, X = X, beta0 = rep(0, ncol(X)),
                            constraint_type = "beta"))
# solve_scaled_lasso(y = y, X = X, beta0 = rep(0, ncol(X)), lambda = 10^-1,
#                    use_constraint = TRUE, constraint_type = "yhat")
# solve_scaled_lasso(y = y, X = X, beta0 = beta, lambda = 10^-1,
#                    use_constraint = TRUE, constraint_type = "yhat")

## Solve multivariate also gives very similar results when using the
## proper lambda (as expected)
# solve_multivariate(beta0 = Beta0, y = y, X = X, lambda = nrow(X) * 10^-1,
#                    constraint_type = "yhat")


test_that("solve_scaled_lasso() has correct dimensions", {
  expect_length(scla, 6)
  expect_length(scla$par$estimate, m)
  expect_length(scla$sigma_scaledlasso, 1)
  expect_named(scla, c("par", "sigma_scaledlasso", "method", "value",
                       "iterations", "last_progress"))
  # expect_true(all(scla$beta_init <= 0))
})

scosys <- calculate_Z(X = X)
# scosys <- score_system(X = X, y = y, beta_init = scla$par$estimate,
#                        hsigma = scla$sigma_scaledlasso)

test_that("score_system() has correct dimensions", {
  expect_equal(dim(scosys), dim(X))
})

tau <- noise_factor(X = X, score_system = scosys)

test_that("noise_factor() has correct dimensions", {
  expect_length(tau, m)
})

beta <- update_beta(X = X, y = y, beta_init = scla$par$estimate,
                    score_system = scosys)

test_that("noise_factor() has correct dimensions", {
  expect_length(beta, m)
})

hci <- size_half_confint_shifts(noise_factor = tau,
                                hsigma = scla$sigma_scaledlasso)$half_size

test_that("size_half_confint_shifts() has correct dimensions", {
  expect_length(hci, m)
})

data.frame(lower = beta - hci, estimate = beta,
           upper = beta + hci, signif = abs(beta) > hci)

# bhat <- lasso.proj(x = X, y = y,
#                    betainit = "scaled lasso",
#                    return.Z = TRUE, standardize = FALSE)$bhat
#
# (bhat - beta)
# (bhat - beta) / beta
# plot(bhat, beta)

V <- covariance_noise_matrix(X, scosys)
V
