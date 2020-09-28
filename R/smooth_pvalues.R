#' Hierarchical correction of p-values
#'
#' @param pvalues Vector of p-values.
#' @param esign Vector of e-signs (optional).
#' @param tree Tree, as a \code{phylo} object.
#' @param arg_shiftpunct Arguments to be passed to
#' \code{estimate_shifts} function.
#' @param arg_shiftconf Arguments to be passed to
#' \code{estimate_confint} function.
#'
#' @return A corrected vector.
#' @export
#'
#' @examples
#' pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
#' tree <- force_ultrametric(chlamydiae$tree)
#' smooth_pvalues(pvalue = pval_obs, tree = tree,
#'                arg_shiftpunct = list(alpha = c(0.1, 2),
#'                                      method = "scaled lasso"))
smooth_pvalues <- function(pvalues, esign = NULL, tree,
                            arg_shiftpunct =
                              list(alphaOU = 1, method = "scaled lasso"),
                            arg_shiftconf =
                              list(alpha_conf = 0.05, method = "scoresystem")){
  zscores <- p2z(pvalues, esign)

  arg_shiftpunct <- c(zscores = list(zscores), tree = list(tree), arg_shiftpunct)

  shiftpunct <- do.call(estimate_shifts, arg_shiftpunct)

  arg_shiftconf <- c(shiftpunct = list(shiftpunct), arg_shiftconf)

  shiftconf <- do.call(estimate_confint, arg_shiftconf)

  pull_pvalues(shiftconf)
}
