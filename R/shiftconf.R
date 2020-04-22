#' 'shiftconf' object
#'
#' @param x a list.
#' @param shiftestim a \code{shiftestim} object.
#'
#' @return a \code{shiftconf} object.
#' @export
#'
as_shiftconf <- function(x, shiftestim){

  required_names <- c("shifts_est", "zscores_est", "alpha_conf", "noise_factor",
                      "covariance_noise_matrix", "method")

  if (!is.list(x) ||
      !all(required_names %in% names(x))) {
    stop("'x' must be the output of a confint function.")
  }

  x <- list(shifts_est = x$shifts_est,
            zscores_est = x$zscores_est,
            alpha_conf = x$alpha_conf,
            method = x$method,
            noise_factor = x$noise_factor,
            covariance_noise_matrix = x$covariance_noise_matrix,
            shiftestim = shiftestim,
            optim_info = x[setdiff(names(x), required_names)])

  class(x) <- "shiftconf"
  return(x)
}

#' @rdname as_shiftconf
#'
#' @name print.shiftconf
#'
#' @param x a 'shiftconf' object.
#' @param ... further arguments to be passed to or from other methods.
#' @inheritParams plot_shifts
#'
#' @importFrom utils head
#' @importFrom stringr str_pad
#' @importFrom ape is.binary multi2di
#'
#' @export
print.shiftconf <- function(x, digits = 3, ...){

  txt_tree1 <- paste0("Tree is", ifelse(x$shiftestim$is_bin, " ", " not "),
                      "binary")
  tree <- x$shiftestim$tree

  txt_tree2 <- paste0("with ", length(tree$tip.label), " leafs and ",
                      nrow(tree$edge), " branches\n")
  cat(txt_tree1, txt_tree2)
  cat(paste("Method:", x$method, "\n"))
  cat(paste("Confidence threshold:", x$alpha_conf, "\n"))
  print(head(x$zscores_est))
}


#' @rdname as_shiftconf
#' @inheritParams plot_shifts
#'
#' @importFrom ggstance geom_pointrangeh
#'
#' @export
plot.shiftconf <- function(x, digits = 3, ...){
  p <- plot_shifts(x$shiftestim$tree, shifts = x$shifts_est$estimate,
                   digits = digits, obs_scores = x$shiftestim$zscores_obs)
  # Ok for the moment, must be completly integrated in plot_shifts() later
  facet_plot(p, panel   = "Estimated confidence interval",
             data    = x$zscores_est,
             geom    = geom_pointrangeh,
             mapping = aes(x = .data$estimate, xmin = .data$lower,
                           xmax = .data$upper, color = .data$estimate >= 0,
                           shape = .data$lower * .data$upper > 0),
             show.legend = FALSE)
}
