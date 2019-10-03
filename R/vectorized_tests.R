#' Vectorized tests
#'
#' @param X matrix with samples in columns and taxa in rows
#' @param Y vector of condition
#'
#' @return A list with a \code{p.value} component and eventually a \code{e.sign} component.
#' @export
#' @importFrom attempt stop_if_not
#' @importFrom stats wilcox.test
#'
#' @examples
test_wilcoxon <- function (X, Y) {
  Y <- factor(Y)
  if(length(levels(Y)) != 2){
    stop("You need exactly two different conditions for a Wilcoxon test.")
  }
  Y <- as.numeric(Y)
  obj <- apply(X, 1, function(x) {
    p.value <- suppressWarnings(wilcox.test(x ~ Y)$p.value)
    e.sign <- sign(mean(x[Y == 2]) - mean(x[Y == 1]))
    c(p.value, e.sign)
  })
  return(list(p.value = obj[1, ], e.sign = obj[2, ]))
}

#' @rdname test_wilcoxon
#' @export
#' @importFrom stats kruskal.test
test_kruskalwallis <- function (X, Y) {
  Y <- as.numeric(factor(Y))
  obj <- apply(X, 1, function(x) {
    suppressWarnings(kruskal.test(x ~ Y)$p.value)
  })
  return(list(p.value = obj))
}

#' @rdname test_wilcoxon
#' @export
#' @importFrom stats lm anova
test_fisher <- function (X, Y) {
  Y <- factor(Y)
  obj <- apply(X, 1, function(x) {
    model <- suppressWarnings(lm(x ~ Y))
    anova(model)$`Pr(>F)`[1]
  })
  return(list(p.value = obj))
}
