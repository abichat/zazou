#' Update confidence interval
#'
#' @param x A shiftconf object.
#' @inheritParams as_shiftconf
#'
#' @return a 'shiftestim' object
#' @export
#'
update_confint <- function(x, alpha_conf){

  if(!inherits(x, "shiftconf")){
    stop("x must be a 'shiftconf' object.")
  }

  add_ci_pv(x, alpha_conf = alpha_conf)
}

add_ci_pv <- function(x, alpha_conf){

  if(!inherits(x, "shiftconf")){
    stop("x must be a 'shiftconf' object.")
  }

  noise_mat <- x$covariance_noise_matrix
  noise_fact <- sqrt(diag(noise_mat))
  hsigma <- x$shiftestim$optim_info$sigma_scaledlasso
  mat_incidence <- incidence_matrix(x$shiftestim$tree)

  shifts_est <- df_confint_pvalue(estimate = x$shifts_est$estimate,
                                  sigma = hsigma, tau = noise_fact,
                                  alpha_conf = alpha_conf)

  zscores_est <- df_conf_leaves(shifts = x$shifts_est$estimate,
                                covariance_noise_mat = noise_mat,
                                mat_incidence = mat_incidence,
                                sigma = hsigma, alpha_conf = alpha_conf)

  x$shifts_est <- shifts_est
  x$zscores_est <- zscores_est
  x$alpha_conf <- alpha_conf

  return(x)
}

#' Confidence interval and p-values for estimates
#'
#' From a vector of estimates computes the lower and upper bound and a p-value.
#'
#' The confidence interval is bivariate whereas the p-value is univariate.
#'
#' @param estimate Shift estimate from \code{shiftestim} object.
#' Eventually named.
#' @param sigma Associated standard error, length \code{1}.
#' @param tau Noise factor, same length as \code{estimate}.
#' @param alpha_conf Confidence level.
#'
#' @importFrom stats pnorm qnorm
#' @return A dataframe with 4 (or 5 if \code{estimate} has names) columns.
#' @export
#'
df_confint_pvalue <- function(estimate, sigma, tau, alpha_conf = 0.05){

  names <- names(estimate)
  if(is.null(names)){
    df <- data.frame(estimate = estimate)
  } else {
    estimate <- unname(estimate)
    df <- data.frame(name = names, estimate = estimate,
                     stringsAsFactors = FALSE)
  }

  sigma_tau <- sigma * tau
  half_confint_size <- qnorm(1 - alpha_conf / 2) * sigma_tau # bivariate
  p_value <- pnorm(estimate / sigma_tau) # univariate


  df$lower <- estimate - half_confint_size
  df$upper <- estimate + half_confint_size
  df$pvalue <- p_value

  return(df)

}

#' @rdname df_confint_pvalue
#'
#' @param shifts Punctual shift estimation from desparsified procedure.
#' @param covariance_noise_mat Covariance noise matrix.
#' @param mat_incidence Incidence matrix
#'
#' @export
#'
df_conf_leaves <- function(shifts, covariance_noise_mat, mat_incidence,
                           sigma, alpha_conf = 0.05){

  n <- nrow(mat_incidence)

  tau <- rep(NA, n)

  for(i in seq_len(n)){
    a <- mat_incidence[i, ]
    tau[i] <- sum(crossprod(a, covariance_noise_mat) * a)
  }
  tau <- sqrt(tau)

  zscores_est <- as.numeric(mat_incidence %*% shifts)
  names(zscores_est) <- rownames(mat_incidence)

  df <- df_confint_pvalue(zscores_est, sigma = sigma, tau = tau,
                          alpha_conf = alpha_conf)
  colnames(df)[1] <- "leaf"
  df
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

#' @rdname update_beta_scoresystem
#' @param colwiseinverse The columnwise inverse.
#'
#' @export
update_beta_colwiseinverse <- function(X, y, beta_init, colwiseinverse){
  res <- y - X %*% beta_init
  correction <- colwiseinverse %*% t(X) %*% res / nrow(X)
  as.numeric(beta_init + correction)
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

#' @rdname noise_factor_scoresystem
#' @inheritParams update_beta_colwiseinverse
#' @param XTXn Matrix.
#'
#' @export
noise_factor_colwiseinverse <- function(X, colwiseinverse, XTXn) {
  MXTXxM <- colwiseinverse %*% XTXn %*% t(colwiseinverse)
  sqrt(diag(MXTXxM) / nrow(X))
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


#' @rdname covariance_noise_matrix_scoresystem
#' @inheritParams noise_factor_colwiseinverse
#' @export
#'
covariance_noise_matrix_colwiseinverse <- function(X, colwiseinverse, XTXn){
  colwiseinverse %*% XTXn %*% t(colwiseinverse) / nrow(X)
}
