#' Node to tips maximum distance.
#'
#' Identify the maximum distance from an internal node amongst all
#' tips.
#' @param subtree defined by main script
#' @param igraph.tree a tree converted into graph format by igraph
#' @param node that is the root of the subtree
#' @param tree Newick formatted ACCTRANS tree
#' @param ntips Number of isolates on the phylogenetic tree
#' @import ape

node_tips_max_dist <- function(node, subtree, igraph.tree, tree, ntips){

  # Retrieve tips under said node
  tips_under_node <- subtree[which(subtree <= ntips)]

  # Create empty vector for  storage
  tips_distance_node <- rep(0, length(tips_under_node))
  interation_dist_vector <- 0

  # Distances of tips from node
  for (x in tips_under_node){
    tip <- as.numeric(as.character(x))
    node_tip_dist <- get_path_distance(igraph.tree,
                                       node,
                                       tip,
                                       tree)

    # Going to store the tip distances into a vector
    interation_dist_vector <- interation_dist_vector + 1

    tips_distance_node[interation_dist_vector] <- node_tip_dist

  }

  # Return the max distance
  return(max(tips_distance_node))
}
