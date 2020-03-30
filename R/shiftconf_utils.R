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

  if(x$method == "desparsified"){
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

#' Extract significant indices
#'
#' Extract indices where lower and upper components of the dataframe
#' are both positive or both negative.
#'
#' @param cidf dataframe with \code{"estimate"}, \code{"lower"},
#' \code{"upper"} columns.
#' @inheritParams extract_significant_leaves
#'
extract_from_cidf <- function(cidf, side = c("left", "both", "right")) {

  stopifnot(c("lower", "upper") %in% colnames(cidf))
  stopifnot(cidf$lower <= cidf$upper)
  side <- match.arg(side)

  signif_left <- c()
  signif_right <- c()

  if (side %in% c("left", "both")) {
    signif_left <- which(cidf$upper < 0)
  }
  if (side %in% c("both", "right")) {
    signif_right <- which(cidf$lower > 0)
  }

  return(sort(union(signif_left, signif_right))
  )
}


#' Extract significant leafs
#'
#' @inheritParams as_shiftestim
#' @param side The side where z-scores are significant (\code{"left"}
#' (default), \code{"both"} or \code{"right"}). Used when method is
#' \code{"desparsified lasso"}.
#'
#' @return The names of the significant leafs
#' @export
#'
extract_significant_leaves <- function(x, side = c("left", "both", "right")){
  if(!inherits(x, "shiftestim")){
    stop("x must be a 'shiftestim' object.")
  }
  if(x$method == "desparsified lasso"){

    signif_ind <- extract_from_cidf(x$zscores_est, side = side)
    return(x$zscores_est$leaf[signif_ind])

  } else {
    msg <- paste0("There is no extraction for this method (",
                  x$method, ").")
    stop(msg)
  }
  return(x)
}
