#' Estimate threshold of rejection for t-scores
#'
#' @param t_scores Vector of t-scores.
#' @param alpha Threshold of detection.
#'
#' @return The threshold of rejection for t-scores.
#' @export
#'
estimate_tstar <- function(t_scores, alpha) {
  p <- length(t_scores)
  t_max <- sqrt(2 * log(p) - 2 * log(log(p)))
  t_scores <- sort(t_scores[t_scores <= t_max])
  values <- p * (1 - pnorm(t_scores)) / seq_along(t_scores)
  if (all(values > alpha)) {
    message("t_star is not feasible, falling back to default value.")
    return(sqrt(2 * log(p)))
  }
  min(t_scores[values <= alpha])
}

#' Correct p-values
#'
#' @details Correction is applicable at a fixed threshold.
#'
#' @param pv Vector of p-values to correct
#' @param alpha Threshold of detection
#'
#' @return The corrected p-values
#' @export
#'
correct_pvalues <- function(pv, alpha) {
  ts <- qnorm(pv)
  t_star <- estimate_tstar(ts, alpha = alpha)
  qv <- pmin(1, pv * 0.05 / pnorm(-t_star))
  names(qv) <- names(pv)
  return(qv)
}
