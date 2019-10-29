
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
  N <- sum(est_shifts != 0)
  - 2 * likelihood(obs_zscores, est_zscores, sigma) + log(N) / 2
}

#' @rdname bic
#' @importFrom stats dnorm
likelihood <- function(obs_zscores, est_zscores, sigma){
  prod(dnorm(x = obs_zscores, mean = est_zscores, sd = sigma))
}
