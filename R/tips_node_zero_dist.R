#' tips_node_zero_dist
#'
#' Identify tips within zero distance of internal node
#' @param subtree Subtree that is made by the main script
#' @param igraph.tree a tree converted into graph format by igraph
#' @param node that is the root of the subtree
#' @param tree Newick formatted ACCTRANS tree
#' @param ntips Number of isolates on the phylogenetic tree
#' tips_node_zero_dist()
tips_node_zero_dist <- function(node, subtree, igraph.tree, tree, ntips){
  tips_under_node <- subtree[which(subtree <= ntips)] #retrieve tips under said node
  tips_zero_distance <- rep(0, length(tips_under_node)) #Create empty vector for storage
  interation_0dist_vector <- 0

  #distances of tips from node
  for (x in tips_under_node){
    tip <- as.numeric(as.character(x))
    node_tip_dist <- get_path_distance(igraph.tree, #graph
                                       node, #node that is the root of the subtree
                                       tip, #tip define distance from said node
                                       tree) #retrieving the distances from node to tip
    #going to store the tips that have a zero into a vector
    if (node_tip_dist == 0){
      interation_0dist_vector <- interation_0dist_vector + 1
      tips_zero_distance[interation_0dist_vector] <- tip
    }
  }
  tips_zero_distance <- tips_zero_distance[tips_zero_distance != 0] #clearing elements in the vector that are zero
  return(tips_zero_distance)
}
