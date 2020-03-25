as_shiftconf <- function(x, shiftestim){

  required_names <- c("shifts_est", "zscores_est", "alpha_conf", "method")

  if (!is.list(x) ||
      !all(required_names %in% names(x))) {
    stop("'x' must be the output of a confint function.")
  }

  x <- list(shifts_est = x$shifts_est,
            zscores_est = x$zscores_est,
            alpha_conf = x$alpha_conf,
            method = x$method,
            shiftestim = shiftestim,
            optim_info = x[setdiff(names(x), required_names)])

  class(x) <- "shiftconf"
  return(x)
}
