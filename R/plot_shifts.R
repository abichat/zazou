#' Plot shift amplitude and location on the tree
#'
#' @param shifts Branch-associated shifts vector, sorted in cladewise order
#' @param tree   Phylo-class object
#' @param scores Leaf-associated scores, as a named vector. If scores is unnamed, assumes scores are in
#'               the same order as \code{tree$tip.label}
#'
#' @return a ggplot object, as created by ggtree
#' @export
#'
#' @importFrom dplyr filter
#' @import ggtree
#'
#' @examples
#' tree <- read.tree(text = "(((A,B),C),(D,E));")
#' shifts <- c(0, -2.9513814680442, 0, 0, 0, -0.00195465052163029, -1.896800199013, 0)
#' scores <- c(-3, -3, 0, -2, 0)
#' plot_shifts(shifts, tree, scores)
plot_shifts <- function(shifts, tree, scores = NULL) {
  ## reorder tree
  tree <- reorder(tree, order = "cladewise")
  ## Build and add branch annotation
  edge_data <- data.frame(node        = tree$edge[, 2],
                          shift_value = shifts) %>%
    dplyr::filter(shift_value != 0)
  p <- ggtree(tree) %<+% edge_data +
    geom_tiplab(size = 5) +
    geom_label(aes(x = branch, label = shift_value))
  ## Build and add zscores annotation
  if (!is.null(scores)) {
    if (is.null(names(scores))) names(scores) <- tree$tip.label
    leaf_data <- data.frame(label = names(scores),
                            score = unname(scores))
    p <- facet_plot(p,
                    panel   = "Z-scores",
                    data    = leaf_data,
                    geom    = geom_point,
                    mapping = aes(x = score))
  }
  p
}
