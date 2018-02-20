#' maria
#'
#' Grouping isolates together using an ACCTRANS tree, from root-to-tip.
#' @param tree Newick formatted ACCTRANS tree.
#' @param thresh SNP threshold for sub-grouping isolates.
#' @param rthreshold The threhold number of ancestral nodes to the root
#' each sub-group must have if they are to be declared as major sub-groups
#' @import geiger
#' @import ape
#' @import igraph
#' @import BMhyd
#' @import phangorn
#' @import phytools


maria <- function(tree,thresh,rthreshold){

  # Retrieving date for output file naming purposes.
  Date <- Sys.Date()

  tree <- read.newick(tree)

  # Resolving dichotomies into a multichotomies.
  tree <- di2multi(tree, 0.5)
  tree_name <- paste("NewickTreefile_", Date, ".tre")
  write.tree(tree, file = tree_name) # Output the new tree.

  # Display Number of Taxa in Tree.
  cat("\t", length(tree$tip.label), " taxa found from phylogenetic tree\n", sep = "")

  # Renaming the node labels.
  tree$node.label <- paste("Node", seq(1, tree$Nnode, 1), sep = "_")

  #===========================================================================#
  #                                                                           #
  #                  Define variables for Sub-grouping                        #
  #                                                                           #
  #===========================================================================#

  # Number of tips
  ntips <- Ntip(tree)

  # Convert tree in igraph form - Takes a 2 column matrix edge list
  igraph.tree <- graph.edgelist(tree$edge)


  # Perform Depth First Search on graph
  dfs <- graph.dfs(igraph.tree,
                   root = ntips + 1,
                   neimode = "out",
                   order = TRUE,
                   dist = TRUE)

  #===========================================================================#
  #                                                                           #
  #                       Investigating edge distances                        #
  #                                                                           #
  #===========================================================================#

  assign <- subgrp(dfs, tree, assign, igraph.tree, thresh)

  sgnum <- max(assign)

  #===========================================================================#
  #                                                                           #
  #                   Identification of Major Sub-Grouping                    #
  #                                                                           #
  #===========================================================================#

  majSubs <- majorsub(sgnum, tree, assign, rthreshold)

  numbmajorsubgroup <- max(majSubs$maj_subgroup_assign)

  #===========================================================================#
  #                                                                           #
  #                   Collating analysis for output                           #
  #                                                                           #
  #===========================================================================#

  # Number of singletons after sub-grouping
  remaining_singletons <- which(assign[1:ntips] == 0)

  data <- collectdata(tree,assign,majSubs)

   # collating variables from above
   output <- list(

     # 1: Number of Subgroups Identified
     SGNo = sgnum,

     # 2: Number of Major Subgroups
     majorSBno = numbmajorsubgroup,

     # 3: Number of Singletons remaining
     singletons = remaining_singletons,

     # 4: Output table of analysis
     output = data,

     # 5: Stating the number of tips in the tree
     ntips = ntips
   )

   cat("",
       paste("Number of Isolates on tree: ", ntips),
       paste("Number of Sub-groups identified: ", sgnum),
       paste("Number of Major Subgroups identified: ", numbmajorsubgroup),
       paste("Number of Singletons remain: ", length(remaining_singletons)),
       sep = "\n")

   #Display Sub-Grouping stats to terminal
   for (i in 1:numbmajorsubgroup){
     numbmajorsubgroupmems <- length(which(data[, 3] == i))
     numbmajorsubgroupsgnum <- length(which(majSubs$majorsublist[, 2] == i))
     cat(paste("Major Sub-Group ", i,
               "is composed of ", numbmajorsubgroupsgnum,
               " Sub-Groups & ", numbmajorsubgroupmems,
               " isolates."), sep = "\n")
   }

   return(output)

}
