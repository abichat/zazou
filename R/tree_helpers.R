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
#' @param phy tree
#' @param states named state vector (typically the estimated z-scores)
#'
#' @return Fitch parsimony score
#' @export
#'
#' @examples
#' nwk <- "(((t4:8,(t6:5,t5:5):3):66,t2:74):40,(t7:39,(t1:10,t3:10):29):75);"
#' tree <- ape::read.tree(text = nwk)
#' states <- c(t1 = "A", t2 = "B", t3 = "A", t4 = "C",
#'             t5 = "C", t6 = "C", t7 = "A")
#' fitch(tree, states) ## 2
fitch <- function(phy, states) {
  ## Check state vector
  if (is.null(names(states))) {
    warning("State vector is unnamed, assuming same order as tip labels")
    names(states) <- phy$tip.label
  }
  ## Reorder tree for up-recursion and state vector to match tip_label order
  phy <- reorder(phy, order = "postorder")
  states <- states[phy$tip.label]
  ## Use integer coding for states for faster subsets
  ## \code{states[tip]} is the integer coding of the state of \code{tip}
  states <- as.integer(factor(states))
  ## Bookkeeping variables
  n_leaves <- length(phy$tip.label)
  n_nodes <- phy$Nnode + n_leaves
  n_states <- length(unique(states))
  costs <- vector(mode = "list", length = n_nodes)
  ## Local functions
  init_cost <- function(node) {
    ## Initialize costs at 0 for internal nodes
    if (node > n_leaves) {
      cost <- rep(0, n_states)
    } else {
      ## Initialize costs at +Inf/0 for leaves according to leaf states
      cost <- rep(Inf, n_states)
      cost[states[node]] <- 0
    }
    costs[[node]] <<- cost
  }
  exists_cost <- function(node) { !is.null(costs[[node]]) }
  update_cost <- function(parent, child) {
    if (!exists_cost(parent)) init_cost(parent)
    if (!exists_cost(child)) init_cost(child)
    costs[[parent]] <<- costs[[parent]] + pmin(
      costs[[child]],         ## parent in same state as child
      min(costs[[child]]) + 1 ## parent in different state
    )
    ## remove child costs
    costs[[child]] <<- numeric(0)
  }
  ## Recursion
  for (edge in seq_len(nrow(phy$edge))) {
    update_cost(parent = phy$edge[edge, 1], child = phy$edge[edge, 2])
  }
  ## Compute fitch scores (min of root score)
  return(min(costs[[n_leaves + 1]]))
}
