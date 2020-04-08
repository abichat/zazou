
#' Information criteria
#'
#' Loglikelihood, Bayesian Information Criterion and Phylogenetic BIC
#'
#' @param obs_zscores Observed z-scores.
#' @param est_zscores Estimated z-scores.
#' @param est_shifts Estimated z-scores.
#' @param mat_covar Covariance matrix.
#'
#' @return scalar
#' @export
#'
#' @examples
#' bic(obs_zscores = c(-2.7, -2.3, 0.2, -0.1, -2, 0.1),
#'     est_zscores = c(-2.5, -2.5, 0, 0, -1, 0),
#'     est_shifts = c(0, 0, 0, 0, -1, 0, -2.5, 0, 0, 0),
#'     mat_covar = diag(rep(1, 6)))
bic <- function(obs_zscores, est_zscores, est_shifts, mat_covar){
  LL <- loglikelihood(obs_zscores = obs_zscores,
                      est_zscores = est_zscores,
                      mat_covar = mat_covar)
  k <- norm0(est_shifts)
  N <- length(obs_zscores)
  - 2 * LL + k * log(N)
}

#' @rdname bic
#' @importFrom mnormt dmnorm
loglikelihood <- function(obs_zscores, est_zscores, mat_covar){
  dmnorm(x = obs_zscores, mean = est_zscores, varcov = mat_covar, log = TRUE)
}

#' @rdname bic
#' @param mat_incidence Incidence matrix.
#' @export
pbic <- function(obs_zscores, est_zscores, est_shifts,
                 mat_incidence, mat_covar){
  LL <- loglikelihood(obs_zscores, est_zscores, mat_covar = mat_covar)
  pos_shifts <- est_shifts != 0
  k <- sum(pos_shifts)
  N <- length(obs_zscores)
  nu <- as.vector(var(obs_zscores))
  res_inc <- mat_incidence[, pos_shifts]
  inv_cov <- solve(mat_covar)
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
#' @param ... Not used, here for avoiding errors when passing unknown arguments
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
