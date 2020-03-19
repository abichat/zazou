#' Update confidence interval
#'
#' @param x a 'shiftestim' object.
#' @param alpha_confint Confidence level.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return a 'shiftestim' object
#' @export
#'
update_confint <- function(x, alpha_confint, ...){

  if(!inherits(x, "shiftestim")){
    stop("x must be a 'shiftestim' object.")
  }

  if(x$method == "desparsified lasso"){
    tau <- x$optim_info$noise_factor
    V <- x$optim_info$covariance_noise_matrix
    hsigma <- x$optim_info$hsigma_scaledlasso
    mat_incidence <- incidence_matrix(x$tree)

    shcs <- size_half_confint_shifts(noise_factor = tau,
                                     hsigma = hsigma,
                                     alpha_conf = alpha_confint)$half_size
    shcz <- size_half_confint_zscores(covariance_noise_mat = V,
                                      incidence_mat = mat_incidence,
                                      hsigma = hsigma,
                                      alpha_conf = alpha_confint)

    x$optim_info$alpha_confint <- alpha_confint
    x$shift_est$lower <- x$shift_est$estimate - shcs
    x$shift_est$upper <- x$shift_est$estimate + shcs
    x$zscores_est$lower <- x$zscores_est$estimate - shcz
    x$zscores_est$upper <- x$zscores_est$estimate + shcz

  } else {
    msg <- paste0("There is no confindence interval for this method (",
                  x$method, ").")
    warning(msg)
  }
  return(x)
}
