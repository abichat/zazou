#' Update confidence interval
#'
#' @inheritParams as_shiftestim
#' @param alpha_confint Confidence level.
#'
#' @return a 'shiftestim' object
#' @export
#'
update_confint <- function(x, alpha_confint){

  if(!inherits(x, "shiftconf")){
    stop("x must be a 'shiftconf' object.")
  }

  if(x$method == "desparsified lasso"){
    tau <- x$optim_info$noise_factor
    V <- x$optim_info$covariance_noise_matrix
    hsigma <- x$shiftestim$optim_info$sigma_scaledlasso
    mat_incidence <- incidence_matrix(x$shiftestim$tree)

    shcs <- size_half_confint_shifts(noise_factor = tau,
                                     hsigma = hsigma,
                                     alpha_conf = alpha_confint)$half_size
    shcz <- size_half_confint_zscores(covariance_noise_mat = V,
                                      incidence_mat = mat_incidence,
                                      hsigma = hsigma,
                                      alpha_conf = alpha_confint)

    x$alpha_conf <- alpha_confint
    x$shifts_est$lower <- x$shifts_est$estimate - shcs
    x$shifts_est$upper <- x$shifts_est$estimate + shcs
    x$zscores_est$lower <- x$zscores_est$estimate - shcz
    x$zscores_est$upper <- x$zscores_est$estimate + shcz

  } else {
    msg <- paste0("There is no confindence interval for this method (",
                  x$method, ").")
    warning(msg)
  }
  return(x)
}


#' Confidence interval and p-values for estimates
#'
#' @param estimate Punctual estimate. Eventually named
#' @param sigma Associated standard error, length \code{1}.
#' @param tau Noise factor, same length as \code{estimate}.
#' @param alpha_conf Confidence level.
#'
#' @return A dataframe with 4 (or 5 if \code{estimate} has names) columns.
#' @export
#'
df_confint_pvalue <- function(estimate, sigma, tau, alpha_conf){

  names <- names(estimate)
  if(is.null(names)){
    df <- data.frame(estimate = estimate)
  } else {
    estimate <- unname(estimate)
    df <- data.frame(name = names, estimate = estimate)
  }

  sigma_tau <- sigma * tau
  half_confint_size <- qnorm(1 - alpha_conf / 2) * sigma_tau
  p_value <- 2 * (1 - pnorm(estimate / sigma_tau)) # pnorm = distribution/repartition

  df$lower <- estimate - half_confint_size
  df$upper <- estimate + half_confint_size
  df$pvalue <- p_value

  return(df)

}
