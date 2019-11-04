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


#' Force a tree to be ultrametric
#'
#' Force a tree to be ultrametric by extending its terminal branches
#'
#' @param tree a phylo object
#'
#' @return an ultrametric tree
#' @export
#' @importFrom ape vcv Ntip
#'
#' @examples
#' tree_u <- force_ultrametric(alcohol$tree)
#' ape::is.ultrametric(tree_u)
force_ultrametric <- function(tree) {
  h <- diag(vcv(tree))
  d <- max(h) - h
  ii <- vapply(1:Ntip(tree),
               FUN = function(x, y) {which(y == x)},
               FUN.VALUE = integer(1),
               y = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + d
  tree
}


#' Fitch parsimony score
#'
#' @param tree tree
#' @param zscores named z-scores
#'
#' @return Fitch parsimony score
#' @export
#'
#' @importFrom phangorn as.phyDat fitch
parsimony_score <- function(tree, zscores){
  z <- as.phyDat(as.factor(zscores))
  suppressWarnings(fitch(tree, z))
}
