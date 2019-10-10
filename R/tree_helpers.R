#' Tree height
#'
#' Compute the tree height from an ultrametric tree
#'
#' @param tree an ultrametric tree \code{phylo} object.
#'
#' @return the height of the tree
#' @export
#' @importFrom ape is.ultrametric
#'
#' @examples
#' tree <- ape::rcoal(6)
#' tree_height(tree)
tree_height <- function(tree){

  stopifnot(is.ultrametric(tree))

    N <- length(tree$tip.label)

  e <- which(tree$edge[, 2] == 1)
  parent <- tree$edge[e, 1]
  sum <- tree$edge.length[e]

  while (parent > N+1) {
    e <- which(tree$edge[, 2] == parent)
    parent <- tree$edge[e, 1]
    sum <- sum + tree$edge.length[e]
  }
  return(sum)
}
