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


#' Size of the half confidence interval for the shifts
#'
#' Compute the size of the half confidence interval for the shifts,
#' depending on the desired level of confidence \code{alpha}.
#'
#' @param noise_factor Noise factor, size m.
#' @param hsigma Estimate value of sigma, found with scaled lasso.
#' @param alpha_conf The confidence level.
#'
#' @importFrom stats qnorm
#'
#' @return The half-size of the confidence interval for each beta, size m+n.
#' @export
size_half_confint_shifts <- function(noise_factor, hsigma, alpha_conf = 0.05){
  stopifnot(length(hsigma) == 1)
  stopifnot(length(alpha_conf) == 1)

  return(list(half_size = qnorm(1 - alpha_conf / 2) * noise_factor * hsigma,
              alpha_conf = alpha_conf))
}



#' Size of the half confidence interval for the z-scores
#'
#' Compute the size of the half confidence interval for the z-scores,
#' depending on the desired level of confidence \code{alpha}.
#'
#' @param covariance_noise_mat The covariance of the noise component
#' for the leaves.
#' @param incidence_mat The incidence matrix of the tree, size \code{n*(n+m)}.
#' @inheritParams size_half_confint_shifts
#'
#' @return The half-size of the confidence interval for each z-score, size n.
#' @export
#'
size_half_confint_zscores <-
  function(covariance_noise_mat, incidence_mat, hsigma, alpha_conf = 0.05){
    n <- nrow(incidence_mat)
    first_term <- size_half_confint_shifts(noise_factor = rep(1, n),
                                           hsigma = hsigma,
                                           alpha_conf = alpha_conf)$half_size
    second_term <- rep(NA, n)
    for(i in seq_len(n)){
      a <- incidence_mat[i, ]
      second_term[i] <- sum(crossprod(a, covariance_noise_mat) * a)
    }
    second_term <- sqrt(second_term)

    first_term * second_term
  }



