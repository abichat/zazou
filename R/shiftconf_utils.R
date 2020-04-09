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
