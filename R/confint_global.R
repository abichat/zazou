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
