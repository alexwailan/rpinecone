#' Sub-lineage isolates on a phyologentic tree.
#'
#' Function to assign sub-lineages to isolates according to the
#' structure of a rooted tree.
#'
#'@param dfs The depth-first search performed on the igraph tree.
#'@param tree The rooted phylogenetic tree.
#'@param assign Vector of Sub-lineage numbers corresponding to tips on a tree.
#'@param ntips Number of tips on the phylogenetic tree.
#'@param igraph.tree Phylogenetic tree converted into a graph
#'@param thresh SNP threshold for sub-lineage isolates on the phylogentic tree.



sublineage <- function(dfs, tree, assign, igraph.tree, thresh, quiet){

  slnum <- 0
  ntips <- Ntip(tree)

    # generating a vector the size of ntips+tree$Nnode corresponding to the number of positions within the tree data
  assign <- rep(0, ntips + tree$Nnode)

  # Traverse the tree in depth first order - starting at the root
  for (i in 1:length(dfs$order)){
    node <- dfs$order[i]

    # Skip leaves/tips
    if(!quiet){
      if (node < ntips + 1) {
        cat(paste(tree$tip.label[node], " skipped as this is a leaf."), sep = "\n")
      } else {
        #if the number of the "node" is equal or more than ntips+1 it is actually a node; proceed with below
        cat("", paste(tree$node.label[node - ntips], " will be investigated."), sep = "\n")
      }
    }
    if (node < ntips + 1) next

    #Node must no be assigned to a Sub-lineage
    if (assign[node] <= 0){

      #Need to generate a subree for sub-lineage analysis
      subtree <- graph.dfs(igraph.tree, node, neimode = "out", unreachable = FALSE)$order

      subtree <- subtree[!is.na(subtree)]

      #Define maximum distance amongst all tips from said node
      node_tips_max_dist <- node_tips_max_dist(node, subtree, igraph.tree, tree, ntips)

      #Define the tips connected to node at zero distance
      tips_0_dist_node <- tips_node_zero_dist(node, subtree, igraph.tree, tree, ntips)

      #Logical sub-lineage numbering; composed of two parameters - threshold && zero distance
      if (node_tips_max_dist <= thresh){

        slnum <- slnum + 1

        assign[subtree] <- slnum

        if(!quiet){
          cat(paste("Threhold met - Sub-lineage Number assigned: ", slnum), sep = "\n")
        }

      } else
        if (sum(tips_0_dist_node) > 0 &&
            tips_0_dist_node[1] != sum(tips_0_dist_node)){

          # Sub-lineage tips in zero distance to node when max distance is over threshold
          slnum <- slnum + 1

          assign[tips_0_dist_node] <- slnum

          assign[node] <- slnum

          if(!quiet){
            cat(paste("Internal sub-lineage detected - Sub-lineage number assigned: ", slnum), sep = "\n")
          }
        }

    } #End of If statement

  } # End of for loop

  return(assign)

}
