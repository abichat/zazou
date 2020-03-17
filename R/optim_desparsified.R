#' Desparsified lasso
#'
#' @param Delta0 Initial value for Delta.
#' @param Y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param alpha_conf Confidence interval for Delta values.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list composed by \code{par} (beta estimate) among others.
#' @export
#'
solve_desparsified <- function(Delta0, Y, X, alpha_conf = 0.05, ...){
  scla <- solve_scaled_lasso(y = Y, X = X, beta0 = Delta0, ...)

  scosys <- calculate_Z(X = X)

  tau <- noise_factor(X = X, score_system = scosys)

  beta <- update_beta(X = X, y = Y, beta_init = scla$par$estimate,
                      score_system = scosys)

  hci <- size_half_confint_shifts(noise_factor = tau,
                                  hsigma = scla$sigma_scaledlasso, alpha_conf)

  par <- data.frame(estimate = beta,
                    lower = beta - hci$half_size,
                    upper = beta + hci$half_size)

  mat_covar_noise <- covariance_noise_matrix(X = X, score_system = scosys)

  list(par = par, value = NA, method = "desparsified lasso",
       alpha_confint = alpha_conf, covariance_noise_matrix = mat_covar_noise)

}

#' Noise factor
#'
#' @param X Matrix size m*(n+m).
#' @param score_system Score system of X.
#'
#' @return The vector of noise factor for the shifts, size (n+m)
#' @export
noise_factor <- function(X, score_system) {
  num <- sqrt(colSums(score_system ^ 2))
  # den <- colSums(X * score_system)
  den <- nrow(X)
  num / den
}

#' Update beta
#'
#' @param y A vector of size m.
#' @param beta_init Initial value of beta found with scaled lasso.
#' @inheritParams noise_factor
#'
#' @return The one-step self-bias corrected estimator of beta, size m.
#' @export
update_beta <- function(X, y, beta_init, score_system) {
  res <- y - X %*% beta_init
  num <- t(score_system) %*% res
  # den <- colSums(score_system * X)
  den <- nrow(X)
  beta_init + num / den
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
#' @return The half-size of the confidence interval for each beta, size m+n.
#' @export
size_half_confint_shifts <- function(noise_factor, hsigma, alpha_conf = 0.05){
  stopifnot(length(hsigma) == 1)
  stopifnot(length(alpha_conf) == 1)

  return(list(half_size = qnorm(1 - alpha_conf / 2) * noise_factor * hsigma,
              alpha_conf = alpha_conf))
}


#' Covariance noise matrix
#'
#' @inheritParams noise_factor
#'
#' @return The covariance of the noise component for the leaves (size n*n)
#' @export
#'
covariance_noise_matrix <- function(X, score_system){
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
                                           alpha_conf = alpha_conf)
    second_term <- rep(NA, n)
    for(i in seq_len(n)){
      a <- incidence_mat[i, ]
      second_term <- sum(crossprod(a, covariance_noise_mat) * a)
    }
    second_term <- sqrt(second_term)

    first_term * second_term
}


