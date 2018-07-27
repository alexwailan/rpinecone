#' pinecone
#'
#' Sub-lineage identification using an ACCTRANS tree, from root-to-tip.
#'
#' @import geiger
#' @import ape
#' @import igraph
#' @import BMhyd
#' @import phangorn
#' @import phytools
#'
#' @param tree a phylo object
#' @param thresh SNP threshold for sub-lineage identification.
#' @param rthreshold The threhold number of ancestral nodes to the root
#' each sub-lineage must have if they are to be declared as major sub-lineages
#' @param quiet whether or not to print out extra information (default=FALSE)
#'
#' @examples
#' tree.file.name <- system.file("extdata", "pyjar.staph_ST2371_45_re_itol_7079_1_6_root.joint.tre", package = "rPinecone")
#' tree <- phytools::read.newick(tree.file.name)
#' pinecone(tree, 2, 3)
#'
#' @export
pinecone <- function(tree, thresh, rthreshold, quiet=FALSE){

  # Some checks
  if(class(tree)!="phylo") stop("tree must be a phylo object!")
  if(!(is.numeric(thresh) && thresh>0)) stop("thresh must be a positive integer!")
  if(!(is.numeric(rthreshold) && rthreshold>0)) stop("rthreshold must be a positive integer!")

  # Retrieving date for output file naming purposes.
  Date <- Sys.Date()

  # Resolving dichotomies into a multichotomies.
  tree <- di2multi(tree, 0.5)

  if(!quiet){
    # Display Number of Taxa in Tree.
    cat("\t", length(tree$tip.label), " taxa found from phylogenetic tree\n", sep = "")
  }

  # Renaming the node labels.
  tree$node.label <- paste("Node", seq(1, tree$Nnode, 1), sep = "_")

  #===========================================================================#
  #                                                                           #
  #        Define variables for Sub-lineage identification                    #
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

  assign <- sublineage(dfs, tree, assign, igraph.tree, thresh, quiet)

  slnum <- max(assign)

  #===========================================================================#
  #                                                                           #
  #                   Identification of Major Sub-lineages                    #
  #                                                                           #
  #===========================================================================#

  majSubs <- majorlineage(slnum, tree, assign, rthreshold)

  numbmajorsublineage <- max(majSubs$maj_sublineage_assign)

  #===========================================================================#
  #                                                                           #
  #                   Collating analysis for output                           #
  #                                                                           #
  #===========================================================================#

  # Number of singletons after sub-lineage identification
  remaining_singletons <- which(assign[1:ntips] == 0)

  data <- collectdata(tree,assign,majSubs)

   # collating variables from above
   output <- list(

     # 1: Number of Sublineages Identified
     SLno = slnum,

     # 2: Number of Major Sublineages
     majorSLno = numbmajorsublineage,

     # 3: Number of Singletons remaining
     singletons = remaining_singletons,

     # 4: Output table of analysis
     table = data,

     # 5: Stating the number of tips in the tree
     ntips = ntips,

     # 6: The phylogentic tree with dichotomies resolved into multichotomies.
     tree=tree
   )

   if(!quiet){
     cat("",
         paste("Number of Isolates on tree: ", ntips),
         paste("Number of Sub-lineages identified: ", slnum),
         paste("Number of Major Sublineages identified: ", numbmajorsublineage),
         paste("Number of Singletons remain: ", length(remaining_singletons)),
         sep = "\n")

     #Display Sub-lineage identification stats to terminal
     for (i in 1:numbmajorsublineage){
       numbmajorsublineagemems <- length(which(data[, 3] == i))
       numbmajorsublineageslnum <- length(which(majSubs$majorsublist[, 2] == i))
       cat(paste("Major Sub-Lineages ", i,
                 "is composed of ", numbmajorsublineageslnum,
                 " Sub-lineages & ", numbmajorsublineagemems,
                 " isolates."), sep = "\n")
     }
   }

   return(output)

}
