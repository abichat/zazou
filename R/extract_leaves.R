#' Extract significant leafs
#'
#' @param x Object of class
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return The names of the significant leafs.
#' @export
#' @seealso \code{\link{extract_significant_leaves.shiftestim}},
#' \code{\link{extract_significant_leaves.shiftconf}}
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
#' @inheritDotParams extract_significant_leaves
#'
#' @return The names of the significant leafs
#' @export
#'
extract_significant_leaves.shiftestim <-
  function(x, threshold = 0, direction = c("<", "<=", ">", ">="), ...) {

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
#' @param threshold P-value threshold.
#' @inheritDotParams extract_significant_leaves
#'
#' @return The names of the significant leafs.
#' @export
#'
extract_significant_leaves.shiftconf <- function(x, threshold = 0.05, ...) {
    lgl <- x$zscores_est$pvalue < threshold
    x$zscores_est$leaf[lgl]
  }


#' @rdname extract_significant_leaves
#' @export
#'
extract_significant_leaves.NULL <- function(x, ...) {
    NULL
  }

