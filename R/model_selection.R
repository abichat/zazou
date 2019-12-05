
#' Information criteria
#'
#' Loglikelihood, Bayesian Information Criterion and Phylogenetic BIC
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

#' @rdname bic
#' @param alpha alpha from OU.
#' @param tree phylogenetic tree.
#' @export
pbic <- function(obs_zscores, est_zscores, est_shifts, sigma, alpha, tree){
  LL <- loglikelihood(obs_zscores, est_zscores, sigma)
  pos_shifts <- est_shifts != 0
  k <- sum(pos_shifts)
  N <- length(obs_zscores)
  nu <- as.vector(var(obs_zscores))
  inc <- incidence_matrix(tree)
  res_inc <- inc[, pos_shifts]
  cov <- covariance_matrix(tree, alpha)
  inv_cov <- solve(cov)
  mat <- nu * t(res_inc) %*% inv_cov %*% res_inc
  logdet <- determinant(mat, logarithm = TRUE)$modulus[1]
  - 2 * LL + 2 * k * log(2 * N - 3) + 2 * log(N) + logdet
}


#' Helper function to compute a logarithmic lambda grid
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_lambda Length of the lambda grid.
#' @param min_ratio Ratio between the minimum and maximum values in
#' the lambda grid.
#' @param ... Not used, here for avoiding errors when passing unknow arguments
#' from previous ellipsis.
#'
#' @export
#' @return a numeric vector of length \code{nlambda}
#'
lambda_grid <- function(x, y, n_lambda = 10, min_ratio = 1e-2, ...) {
  ## center y and x
  y <- scale(y, center = TRUE, scale = FALSE)
  x <- scale(x, center = TRUE, scale = FALSE)
  xty <- drop(crossprod(x, y))
  lambda_max <- max(abs(xty))
  return(10 ^ seq(from = log10(lambda_max), to = log10(min_ratio * lambda_max),
                  len = n_lambda))
}
