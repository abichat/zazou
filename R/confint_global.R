#' Estimate confidence interval from scaled lasso.
#'
#' @param shiftestim a \code{shiftestim} object with method \code{scaledlasso}.
#' @param alpha_conf the confidence level.
#' @param method For the moment, only \code{scoresystem}.
#' @param ... further argument...
#'
#' @return a \code{shiftconf} object.
#' @export
#'
estimate_confint <- function(shiftestim, alpha_conf = 0.05,
                             method = c("scoresystem", "colwiseinverse"),
                             ...){

  method <- match.arg(method)

  stopifnot(inherits(shiftestim, "shiftestim"))

  obj <- switch (method,

                 "scoresystem" = confint_scoresystem(
                   x = shiftestim, ...),

                 "colwiseinverse" = confint_colwiseinverse(
                   x = shiftestim, ...)

                 )

  obj <- as_shiftconf(obj, shiftestim, alpha_conf = alpha_conf)

  obj$optim_info$supp_arg <- list(...)

  return(obj)
}
