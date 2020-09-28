#' Pull p-values
#'
#' @param x Object of class shiftpunct or shiftconf
#'
#' @return The named vector of p-values
#' @export
#'
pull_pvalues <- function(x){
  UseMethod("pull_pvalues")
}

#' Pull p-values
#'
#' @param x Object of class shiftpunct or shiftconf
#'
#' @return The named vector of p-values
#' @export
#' @importFrom stats pnorm
#'
pull_pvalues.shiftpunct <- function(x){
  pnorm(x$zscores_est)
}

#' Pull p-values
#'
#' @param x Object of class shiftpunct or shiftconf
#'
#' @return The named vector of p-values
#' @export
#'
pull_pvalues.shiftconf <- function(x){
  pv <- x$zscores_est$pvalue
  names(pv) <- x$zscores_est$leaf
  return(pv)
}
