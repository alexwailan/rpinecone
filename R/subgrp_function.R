#' Subgrp
#'
#' Function to assign sub-groups to isolates according to the structure of a rooted tree
#'
#'@param dfs The depth-first search performed on the igraph tree.
#'@param tree The rooted phylogenetic tree.
#'@param assign Vector of Sub-group numbers corresponding to tips on a tree.
#'@param ntips Number of tips on the phylogenetic tree.
#'@param igraph.tree Phylogenetic tree converted into a graph
#'@param thresh SNP threshold for sub-group isolates on the phylogentic tree.
#'
#' subgrp()


subgrp <- function(dfs, tree, assign, igraph.tree, thresh){

  sgnum <- 0
  ntips <- Ntip(tree)

    # generating a vector the size of ntips+tree$Nnode corresponding to the number of positions within the tree data
  assign <- rep(0, ntips + tree$Nnode)

  # Traverse the tree in depth first order - starting at the root
  for (i in 1:length(dfs$order)){
    node <- dfs$order[i]

    # Skip leaves/tips
    if (node < ntips + 1) {
      cat(paste(tree$tip.label[node], " skipped as this is a leaf."), sep = "\n")
      next
    } else {
      #if the number of the "node" is equal or more than ntips+1 it is actually a node; proceed with below
      cat("", paste(tree$node.label[node - ntips], " will be investigated."), sep = "\n")
    }

    #Node must no be assigned to a Sub-Group
    if (assign[node] <= 0){

      #Need to generate a subree for sub-grouping analysis
      subtree <- graph.dfs(igraph.tree, node, neimode = "out", unreachable = FALSE)$order

      subtree <- subtree[!is.na(subtree)]

      #Define maximum distance amongst all tips from said node
      node_tips_max_dist <- node_tips_max_dist(node, subtree, igraph.tree, tree, ntips)

      #Define the tips connected to node at zero distance
      tips_0_dist_node <- tips_node_zero_dist(node, subtree, igraph.tree, tree, ntips)

      #Logical sub-grouping numbering; composed of two parameters - threshold && zero distance
      if (node_tips_max_dist <= thresh){

        sgnum <- sgnum + 1

        assign[subtree] <- sgnum

        cat(paste("Threhold met - Sub-Group Number assigned: ", sgnum), sep = "\n")

      } else
        if (sum(tips_0_dist_node) > 0 &&
            tips_0_dist_node[1] != sum(tips_0_dist_node)){

          # Sub-grouping tips in zero distance to node when max distance is over threshold
          sgnum <- sgnum + 1

          assign[tips_0_dist_node] <- sgnum

          assign[node] <- sgnum

          cat(paste("Internal sub-group detected - Sub-group number assigned: ", sgnum), sep = "\n")

        }

    } #End of If statement

  } # End of for loop

  return(assign)

}
