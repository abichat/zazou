
#' Information criteria
#'
#' Loglikelihood, Bayesian Information Criterion and Phylogenetic BIC
#'
#' @param obs_zscores Observed z-scores.
#' @param est_zscores Estimated z-scores.
#' @param est_shifts Estimated z-scores.
#' @param mat_covarOU Covariance matrix.
#'
#' @return scalar
#' @export
#'
#' @examples
#' bic(obs_zscores = c(-2.7, -2.3, 0.2, -0.1, -2, 0.1),
#'     est_zscores = c(-2.5, -2.5, 0, 0, -1, 0),
#'     est_shifts = c(0, 0, 0, 0, -1, 0, -2.5, 0, 0, 0),
#'     mat_covarOU = diag(rep(1, 6)))
bic <- function(obs_zscores, est_zscores, est_shifts, mat_covarOU){
  LL <- loglikelihood(obs_zscores = obs_zscores,
                      est_zscores = est_zscores,
                      mat_covarOU = mat_covarOU)
  k <- norm0(est_shifts)
  N <- length(obs_zscores)
  - 2 * LL + k * log(N) * log(log(N))
}

#' @rdname bic
#' @importFrom mnormt dmnorm
loglikelihood <- function(obs_zscores, est_zscores, mat_covarOU){
  dmnorm(x = obs_zscores, mean = est_zscores, varcov = mat_covarOU, log = TRUE)
}

#' @rdname bic
#' @param mat_incidence Incidence matrix.
#' @export
pbic <- function(obs_zscores, est_zscores, est_shifts,
                 mat_incidence, mat_covarOU){
  LL <- loglikelihood(obs_zscores, est_zscores, mat_covarOU = mat_covarOU)
  pos_shifts <- est_shifts != 0
  k <- sum(pos_shifts)
  N <- length(obs_zscores)
  nu <- as.vector(var(obs_zscores))
  res_inc <- mat_incidence[, pos_shifts]
  inv_cov <- solve(mat_covarOU)
  mat <- nu * t(res_inc) %*% inv_cov %*% res_inc
  logdet <- determinant(mat, logarithm = TRUE)$modulus[1]
  - 2 * LL + 2 * k * log(2 * N - 3) + 2 * log(N) + logdet
}



