#' Incidence matrix
#'
#' Compute the incidence matrix from a tree
#'
#' @param phylo a tree
#'
#' @return the incidence matrix
#' @export
#'
#' @importFrom stats reorder
#'
#' @examples
#' tree <- ape::rtree(3)
#' incidence_matrix(tree)
incidence_matrix <- function(phylo) {
  ## Reorder and keep track
  phy <- reorder(phylo, order = "postorder")
  cor <- corresponding_edges(edges = seq_len(nrow(phy$edge)),
                             from = phylo, to = phy)
  ## Init and recurence
  mat <- init_incidence_matrix(phy)
  mat <- recursion_up(phy, mat, update_incidence_matrix)
  ## Take for each node its parenting branch
  daughters <- phy$edge[, 2]
  mat <- mat[daughters, ]
  ## Return to original tree
  mat <- t(mat[cor, ])
  rownames(mat) <- phy$tip.label
  return(mat)
}

init_incidence_matrix <- function(phy) {
  ntaxa <- length(phy$tip.label)
  mat <- matrix(NA, nrow = 1 + nrow(phy$edge), ncol = ntaxa)
  for (i in 1:ntaxa) {
    mat[i, ] <- 1:ntaxa == i
  }
  return(mat)
}

update_incidence_matrix <- function(daughtersParams, ...) {
  inc <- colSums(daughtersParams)
  return(inc > 0)
}

recursion_up <- function(phy, params, updateUp, ...) {
  if (attr(phy, "order") != "postorder")
    stop("The tree must be in postorder order")
  ## Tree recursion
  e <- 1
  while (e <= nrow(phy$edge)) {
    edge <- phy$edge[e,]
    parent <- edge[1]
    ii <- which(phy$edge[, 1] == parent)
    daughters <- phy$edge[ii, 2]
    params[parent, ] <-
      updateUp(edgesNbr = ii, daughters = daughters, parent = parent,
               daughtersParams = params[daughters, , drop = FALSE], ...)
    e <- ii[length(ii)] + 1
  }
  return(params)
}


corresponding_edges <- function(edges, from, to) {
  mm <- match(from$edge[, 2], to$edge[, 2])
  newEdges <- mm[edges]
  return(newEdges)
}




