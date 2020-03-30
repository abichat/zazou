confint_desparsified <- function(x, alpha_conf = 0.05, ...){
  stopifnot(inherits(x, "shiftestim"))

  mat_covar <- covariance_matrix(x$tree, x$alpha)
  mat_incidence <- incidence_matrix(x$tree)
  R <- inverse_sqrt(mat_covar)
  Y <- R %*% x$zscores_obs
  X <- R %*% mat_incidence
  hsigma <- x$optim_info$sigma_scaledlasso

  scosys <- calculate_Z(X = X)

  tau <- noise_factor(X = X, score_system = scosys)

  new_beta <- update_beta(X = X, y = Y, beta_init = x$shift_est,
                          score_system = scosys)

  mat_covar_noise <- covariance_noise_matrix(X = X, score_system = scosys)

  hcis <- size_half_confint_shifts(noise_factor = tau,
                                   hsigma = hsigma, alpha_conf = alpha_conf)

  hciz <- size_half_confint_zscores(covariance_noise_mat = mat_covar_noise,
                                    incidence_mat = mat_incidence,
                                    hsigma = hsigma, alpha_conf = alpha_conf)

  shifts_est <- data.frame(estimate = new_beta,
                           lower = new_beta - hcis$half_size,
                           upper = new_beta + hcis$half_size)

  zscores_est <- mat_incidence %*% shifts_est$estimate
  zscores_est <- data.frame(leaf = rownames(zscores_est),
                            estimate = zscores_est[, 1],
                            stringsAsFactors = FALSE)
  zscores_est$lower <- zscores_est$estimate - hciz
  zscores_est$upper <- zscores_est$estimate + hciz
  rownames(zscores_est) <- NULL

  list(shifts_est = shifts_est, zscores_est = zscores_est,
       noise_factor = tau, covariance_noise_matrix = mat_covar_noise,
       alpha_conf = alpha_conf, method = "desparsified lasso")
}

