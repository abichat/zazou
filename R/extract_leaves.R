#' Extract significant leafs
#'
#' @param x Object of class
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return The names of the significant leafs.
#' @export
#' @seealso [extract_significant_leaves.shiftestim()]
#' [extract_significant_leaves.shiftconf()]
#'
extract_significant_leaves <- function (x, ...) {
  UseMethod("extract_significant_leaves")
}

#' Extract significant leafs
#'
#' @param x A "shiftestim" object.
#' @param threshold The threshold where selection begins.
#' @param direction The direction of the selection from the threshold.
#' Default to \code{"<"}.
#'
#' @return The names of the significant leafs
#' @export
#'
extract_significant_leaves.shiftestim <-
  function(x, threshold = 0, direction = c("<", "<=", ">", ">=")) {

    direction <- match.arg(direction)

    lgl <- switch(direction,
      "<"  = x$zscores_est <  threshold,
      "<=" = x$zscores_est <= threshold,
      ">"  = x$zscores_est >  threshold,
      ">=" = x$zscores_est >= threshold,
    )
      names(which(lgl))
  }

#' Extract significant leafs
#'
#' @param x A "shiftconf" object.
#' @param side The side where z-scores are significant (\code{"left"}
#' (default), \code{"both"} or \code{"right"}).
#'
#' @return The names of the significant leafs.
#' @export
#'
extract_significant_leaves.shiftconf <-
  function(x, side = c("left", "both", "right")) {
    signif_ind <- extract_from_cidf(x$zscores_est, side = side)
    x$zscores_est$leaf[signif_ind]
  }

#' Extract significant indices
#'
#' Extract indices where lower and upper components of the dataframe
#' are both positive or both negative.
#'
#' @param cidf confidence interval dataframe with \code{"estimate"},
#' \code{"lower"}, \code{"upper"} columns.
#' @inheritParams extract_significant_leaves.shiftconf
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
