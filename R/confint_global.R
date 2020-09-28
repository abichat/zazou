#' Estimate confidence interval from scaled lasso.
#'
#' @param shiftpunct a \code{shiftpunct} object with method \code{scaledlasso}.
#' @param alpha_conf the confidence level.
#' @param method For the moment, only \code{scoresystem}.
#' @param ... further argument...
#'
#' @return a \code{shiftconf} object.
#' @export
#'
estimate_confint <- function(shiftpunct, alpha_conf = 0.05,
                             method = c("scoresystem", "colwiseinverse"),
                             ...){

  method <- match.arg(method)

  stopifnot(inherits(shiftpunct, "shiftpunct"))

  obj <- switch (method,

                 "scoresystem" = confint_scoresystem(
                   x = shiftpunct, ...),

                 "colwiseinverse" = confint_colwiseinverse(
                   x = shiftpunct, ...)

                 )

  obj <- as_shiftconf(obj, shiftpunct, alpha_conf = alpha_conf)

  obj$optim_info$supp_arg <- list(...)

  return(obj)
}
