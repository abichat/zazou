#' Score system
#'
#' Compute the score system of a
#'
#' @param X A vector of size m*(n+m).
#' @param y A vector of size m.
#' @param beta_init Initial value of beta found with scaled lasso.
#' @param hsigma Estimate value of sigma, found with scaled lasso.
#'
#' @return The matrix of score system, same dimension as X.
#' @export
#'
#' @importFrom hdi lasso.proj
score_system <- function(X, y, beta_init, hsigma) {
  obj <-
    suppressWarnings(suppressMessages(
      lasso.proj(x = X, y = y, betainit = beta_init,
                 sigma = hsigma, return.Z = TRUE)
    ))
  obj$Z
}

#' Noise factor
#'
#' @inheritParams score_system
#' @param score_system Score system of X
#'
#' @return The vector of noise factor, size (n+m)
#' @export
noise_factor <- function(X, score_system) {
  num <- colSums(score_system ^ 2)
  den <- colSums(X * score_system)
  num / den
}

#' Compute beta
#'
#' @inheritParams score_system
#' @inheritParams noise_factor
#'
#' @return The one-step self-bias corrected estimator of beta.
#' @export
beta <- function(X, y, beta_init, score_system) {
  res <- y - X %*% beta_init
  num <- colSums(score_system * res)
  den <- colSums(score_system * X)
  beta_init + num / den
}
