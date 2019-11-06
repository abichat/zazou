
#' Bayesian information criterion
#'
#' @param obs_zscores observed z-scores
#' @param est_zscores estimated z-scores
#' @param est_shifts estimated z-scores
#' @param sigma standard error
#'
#' @return scalar
#' @export
#'
#' @examples
#' bic(obs_zscores = c(-2.7, -2.3, 0.2, -0.1, -2, 0.1),
#'     est_zscores = c(-2.5, -2.5, 0, 0, -1, 0),
#'     est_shifts = c(0, 0, 0, 0, -1, 0, -2.5, 0, 0, 0),
#'     sigma = 1.4)
bic <- function(obs_zscores, est_zscores, est_shifts, sigma){
  LL <- loglikelihood(obs_zscores, est_zscores, sigma)
  k <- sum(est_shifts != 0)
  N <- length(obs_zscores)
  - 2 * LL + k * log(N)
}

#' @rdname bic
#' @importFrom stats dnorm
loglikelihood <- function(obs_zscores, est_zscores, sigma){
  sum(log(dnorm(x = obs_zscores, mean = est_zscores, sd = sigma)))
}


#' Helper function to compute a logarithmic lambda grid
#'
#' @param x Design matrix
#' @param y Response vector
#' @param n_lambda Length of the lambda grid
#' @param min_ratio Ratio between the minimum and maximum values in the lambda grid
#'
#' @return a numeric vector of length \code{nlambda}
#'
lambda_grid <- function(x ,y, n_lambda = 6, min_ratio = 1e-5) {
  ## center y and x
  y <- scale(y, center = TRUE, scale = FALSE)
  x <- scale(x, center = TRUE, scale = FALSE)
  xty <- drop(crossprod(x, y))
  lambda_max <- max(abs(xty))
  return(10^seq(log10(lambda_max), log10(min_ratio*lambda_max), len = n_lambda))
}
