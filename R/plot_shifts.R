#' Plot shift amplitude and location on the tree
#'
#' Plot
#'
#' If scores is unnamed, assumes scores are in the same order as
#' \code{tree$tip.label}
#'
#' @param shifts Branch-associated shifts vector, sorted in cladewise order
#' @param tree   Phylo-class object
#' @param true_scores Leaf-associated scores, as a named vector.
#' @param obs_scores Leaf-associated scores, as a named vector.
#' @param est_scores Leaf-associated scores, as a named vector.
#' @param sup_scores A list of lists for additional scores to plot. Each
#' list must contain at least a 'scores" element and eventually a
#' 'title" and 'color' one.
#' @param digits Number of digits to round to.
#'
#' @return a \code{ggplot} object, as created by \code{ggtree}
#' @export
#'
#' @importFrom ggtree ggtree aes %<+% geom_tiplab geom_point geom_label
#' theme_tree2 facet_plot
#' @importFrom rlang .data
#'
#' @examples
#' tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
#' shifts <- c(0, -3, 0, 0, 0, 0, -2, 0)
#' scores <- c(-3.2, -2.8, 0.1, -2.1, -0.1)
#' plot_shifts(tree, shifts, obs_scores = scores)
plot_shifts <- function(tree, shifts, true_scores = NULL,
                        obs_scores = NULL, est_scores = NULL,
                        sup_scores = NULL, digits = 3) {

  ## reorder tree
  tree <- reorder(tree, order = "cladewise")

  ## Build and add branch annotation
  edge_data <- data.frame(node = tree$edge[, 2], shift_value = shifts)
  edge_data <- edge_data[edge_data$shift_value != 0, ]

  p <-
    ggtree(tree) %<+%
    edge_data +
    geom_tiplab(linesize = 5) +
    geom_label(aes(x = .data$branch,
                   label = round(.data$shift_value, digits))) +
    theme_tree2()
# browser()
  ## Build and add zscores annotation
  if (!is.null(true_scores)) {
    if (is.null(names(true_scores))) names(true_scores) <- tree$tip.label
    leaf_data_t <- data.frame(label = names(true_scores),
                              score = unname(true_scores))
    p <- facet_plot(p,
                    panel   = "True z-scores",
                    data    = leaf_data_t,
                    geom    = geom_point,
                    mapping = aes(x = .data$score, color = .data$score >= 0),
                    show.legend = FALSE)
  }

  if (!is.null(obs_scores)) {
    if (is.null(names(obs_scores))) names(obs_scores) <- tree$tip.label
    leaf_data_o <- data.frame(label = names(obs_scores),
                              score = unname(obs_scores))
    p <- facet_plot(p,
                    panel   = "Observed z-scores",
                    data    = leaf_data_o,
                    geom    = geom_point,
                    mapping = aes(x = .data$score, color = .data$score >= 0),
                    show.legend = FALSE)
  }

  if (!is.null(est_scores)) {
    if (is.null(names(est_scores))) names(est_scores) <- tree$tip.label
    leaf_data_e <- data.frame(label = names(est_scores),
                              score = unname(est_scores))
    p <- facet_plot(p,
                    panel   = "Estimated z-scores",
                    data    = leaf_data_e,
                    geom    = geom_point,
                    mapping = aes(x = .data$score, color = .data$score >= 0),
                    show.legend = FALSE)
  }

  if (!is.null(sup_scores)) {

    if(!is.list(sup_scores) | length(sup_scores) == 0) {
      stop("'sup_scores' must be a non-empty list.")
    }

    if(!all(vapply(sup_scores, is.list, FUN.VALUE = logical(1)))) {
      stop("'sup_scores' must be a list of lists.")
    }

    for(i in seq_along(sup_scores)) {
      l_scores <- sup_scores[[i]]
      scores <- l_scores$scores

      if(is.null(scores)){
        stop("Each list in 'sup_scores' must contain a 'scores' element.")
      }

      if (is.null(names(scores))) names(scores) <- tree$tip.label
      if (is.null(l_scores$title)) l_scores$title <- i
      if (is.null(l_scores$color)) l_scores$color <- NA

      leaf_data_s <- data.frame(label = names(scores),
                                scores = unname(scores),
                                col = l_scores$color)
      p <- facet_plot(p,
                      panel   = l_scores$title,
                      data    = leaf_data_s,
                      geom    = geom_point,
                      mapping = aes(x = .data$scores, color = .data$col),
                      show.legend = FALSE)

    }
  }

  p
}
