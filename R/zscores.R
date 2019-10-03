#' Transform p-values to z-scores
#'
#' @param p.value a vector of p-values
#' @param e.sign 	a vector of signs of the effects, taken on value -1 and 1
#' @param tol a numeric value at which the p-value (both ends) will truncate
#'
#' @return a vector of z-values
#' @export
#' @importFrom stats qnorm
#'
#' @examples
#' pv <- runif(5)
#' p2z(pv)
#' p2z(0.5) ## should be 0
p2z <- function (p.value, e.sign = NULL, tol = 1e-15) {
  p.value[p.value <= tol] <- tol
  p.value[p.value >= 1 - tol] <- 1 - tol

  if(is.null(e.sign)) {
    z <- qnorm(1 - p.value)
  } else {
    e.sign[e.sign == 0] <- sample(c(-1, 1), sum(e.sign == 0), replace = TRUE)
    z1 <- qnorm(p.value / 2)
    z2 <- qnorm(1 - p.value / 2)
    z <- ifelse(e.sign > 0, z2, z1)
  }

  return(z)
}
