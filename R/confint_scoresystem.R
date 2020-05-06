#' Desparsified lasso based on score system
#'
#' @param x a 'shiftestim' object.
#' @param alpha_conf the confidence level.
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return a list
#' @export
#'
confint_scoresystem <- function(x, alpha_conf = 0.05, ...){
  stopifnot(inherits(x, "shiftestim"))

  mat_covarOU <- covarianceOU_matrix(x$tree, x$alpha)
  mat_incidence <- incidence_matrix(x$tree)
  R <- inverse_sqrt(mat_covarOU)
  Y <- R %*% x$zscores_obs
  X <- R %*% mat_incidence
  hsigma <- x$optim_info$sigma_scaledlasso

  scosys <- calculate_Z(X = X)

  # tau <- noise_factor_scoresystem(X = X, score_system = scosys)

  new_beta <- update_beta_scoresystem(X = X, y = Y, beta_init = x$shifts_est,
                                      score_system = scosys)

  mat_covar_noise <-
    covariance_noise_matrix_scoresystem(X = X, score_system = scosys)

  # shifts_est <- df_confint_pvalue(estimate = new_beta, sigma = hsigma,
  #                                 tau = tau, alpha_conf = alpha_conf)


  # zscores_est <- df_conf_leaves(shifts = shifts_est$estimate,
  #                               covariance_noise_mat = mat_covar_noise,
  #                               mat_incidence = mat_incidence, sigma = hsigma,
  #                               alpha_conf = alpha_conf)


  list(shifts_est = data.frame(estimate = new_beta),
       zscores_est = NA, covariance_noise_matrix = mat_covar_noise,
       method = "scoresystem",
       scoresystem = scosys)
}


#' Update beta
#'
#' @param X Matrix size m*(n+m).
#' @param y A vector of size m.
#' @param beta_init Initial value of beta found with scaled lasso.
#' @param score_system Score system of \code{X}.
#'
#' @return The one-step self-bias corrected estimator of beta, size m.
#' @export
update_beta_scoresystem <- function(X, y, beta_init, score_system) {
  res <- y - X %*% beta_init
  num <- t(score_system) %*% res
  # den <- colSums(score_system * X)
  den <- nrow(X)
  as.numeric(beta_init + num / den)
}

#' Noise factor
#'
#' @inheritParams update_beta_scoresystem
#'
#' @return The vector of noise factor for the shifts, size (n+m)
#' @export
noise_factor_scoresystem <- function(X, score_system) {
  num <- sqrt(colSums(score_system ^ 2))
  # den <- colSums(X * score_system)
  den <- nrow(X)
  num / den
}

#' Covariance noise matrix
#'
#' @inheritParams noise_factor_scoresystem
#'
#' @return The covariance of the noise component for the leaves (size n*n)
#' @export
#'
covariance_noise_matrix_scoresystem <- function(X, score_system){
  STS <- crossprod(score_system)
  STX <- crossprod(score_system, X)
  n <- ncol(STS)
  V <- matrix(NA, nrow = n, ncol = n)
  for(i in seq_len(n)){
    for(j in i:n){
      # cat(paste0("i = ", i, ", j =", j), sep = "\n")
      V[i, j] <- V[j, i] <- STS[i, j] / abs(STX[i, i] * STX[j, j])
    }
  }
  V
}
