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
#' The confidence interval is bivariate whereas the p-value and q-values are
#' univariate. q-values are threshold-dependent and need to be recomputed when
#' \code{alpha_conf} changes.
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
  if (is.null(names)) {
    df <- data.frame(estimate = estimate)
  } else {
    estimate <- unname(estimate)
    df <- data.frame(name = names,
                     estimate = estimate,
                     stringsAsFactors = FALSE)
  }

  sigma_tau <- sigma * tau
  half_confint_size <- qnorm(1 - alpha_conf / 2) * sigma_tau # bivariate
  p_value <- pnorm(estimate / sigma_tau) # univariate
  q_value <- correct_pvalues(p_value, # univariate and threshold-dependent
                             alpha = alpha_conf)


  df$lower <- estimate - half_confint_size
  df$upper <- estimate + half_confint_size
  df$pvalue <- p_value
  df$qvalue <- q_value

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

  for (i in seq_len(n)) {
    a <- mat_incidence[i,]
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




