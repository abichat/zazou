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
  scla <- scaled_lasso(y = Y, X = X, beta0 = Delta0, ...)

  scosys <- calculate_Z(X = X)

  tau <- noise_factor(X = X, score_system = scosys)

  beta <- update_beta(X = X, y = Y, beta_init = scla$par$estimate,
                      score_system = scosys)

  hci <- size_half_confint(noise_factor = tau,
                           hsigma = scla$sigma_scaledlasso, alpha_conf)

  par <- data.frame(estimate = beta,
                    lower = beta - hci$half_size,
                    upper = beta + hci$half_size)

  list(par = par, value = NA, method = "desparsified lasso",
       alpha_confint = alpha_conf)

}

#' Noise factor
#'
#' @param X Matrix size m*(n+m).
#' @param score_system Score system of X.
#'
#' @return The vector of noise factor, size (n+m)
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

#' Size of the half confidence interval
#'
#' Compute the size of the half confidence interval, depending on the desired
#' level of confidence \code{alpha}.
#'
#' @param noise_factor Noise factor, size m.
#' @param hsigma Estimate value of sigma, found with scaled lasso.
#' @param alpha_conf The confidence level.
#'
#' @return The half-size of the confidence interval for each beta, size m.
#' @export
size_half_confint <- function(noise_factor, hsigma, alpha_conf = 0.05){
  stopifnot(length(hsigma) == 1)
  stopifnot(length(alpha_conf) == 1)

  return(list(half_size = qnorm(1 - alpha_conf / 2) * noise_factor * hsigma,
              alpha_conf = alpha_conf))
}




