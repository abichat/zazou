#' Create clusters
#'
#' @param tree tree whose tips must be clustered.
#' @param N_clusters number of clusters.
#' @param method method to create clusters among \code{"monophyletic"},
#' \code{"paraphyletic"} and \code{"uniform"}.
#'
#' @return The belonging of each tip to a cluster as a named integer.
#'
#' @importFrom stats var cutree cophenetic
#' @importFrom ape as.hclust.phylo
#' @importFrom cluster pam
#' @export
#'
#' @examples
#' tree <- ape::rtree(10)
#' N <- 3
#' plot(tree)
#' create_clusters(tree, N)
#' create_clusters(tree, N, method = "paraphyletic")
create_clusters <- function(tree, N_clusters,
                            method = c("monophyletic", "paraphyletic",
                                       "uniform")){

  method <- match.arg(method)
  method

  N_tips <- length(tree$tip.label)

  if(N_clusters == N_tips){
    clustering <- seq_len(N_clusters)
    names(clustering) <- tree$tip.label
    return(clustering)
  }

  if(N_clusters > N_tips || N_clusters <= 0){
    stop(paste0("Number of clusters clusters must be between 1 and ",
                 N_tips, "."))
  }

  if(method == "uniform"){
    clustering <- sample(x = seq_len(N_clusters),
                         size = N_tips, replace = TRUE)
    while(length(unique(clustering)) != N_clusters){
      clustering <- sample(x = seq_len(N_clusters),
                           size = N_tips, replace = TRUE)
    }
    names(clustering) <- tree$tip.label
  } else if (method == "monophyletic"){
    tree_hcl <- as.hclust.phylo(force_ultrametric(tree))
    clustering <- cutree(tree_hcl, k = N_clusters)
  } else if (method == "paraphyletic"){
    D <- cophenetic(tree)
    partition <- pam(D, N_clusters)
    clustering <- partition$clustering
  }

  return(clustering)
}

