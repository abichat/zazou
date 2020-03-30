#' Estimate confidence interval from scaled lasso.
#'
#' @param shiftestim a \code{shiftestim} object with method \code{scaledlasso}.
#' @param alpha_conf the confidence level.
#' @param method For the moment, only \code{desparsified}.
#' @param ... further argument...
#'
#' @return a \code{shiftconf} object.
#' @export
#'
estimate_confint <- function(shiftestim, alpha_conf = 0.05,
                             method = c("desparsified"), ...){

  stopifnot(inherits(shiftestim, "shiftestim"))

  obj <- switch (method,

                 "desparsified" = confint_desparsified(
                   x = shiftestim, alpha_conf = alpha_conf, ...)

                 )

  obj <- as_shiftconf(obj, shiftestim)

  return(obj)
}
