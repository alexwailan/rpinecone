#' boot_pinecone
#'
#' A wrapper around the main pinecone function that takes multiple trees as input and infers a co-occupancy matrix representing uncertainty in the pinecone clusters.
#'
#'
#' @param tree_list a list of phylogenies or a 'multiPhylo' object of trees that have branch lengths scaled to SNP distance
#' @param thresh SNP threshold for sub-lineage identification.
#' @param rthreshold The threhold number of ancestral nodes to the root
#' each sub-lineage must have if they are to be declared as major sub-lineages
#' @param quiet whether or not to print out extra information (default=TRUE)
#'
#' @return A matrix representing the co-occupancy of isolates in pinecone clusters
#'
#' @examples
#' tree.file.name <- system.file("extdata", "pyjar.staph_ST2371_45_re_itol_7079_1_6_root.joint.tre", package = "rPinecone")
#' tree <- phytools::read.newick(tree.file.name)
#' boot_pinecone(list(tree, tree), 2, 3)
#'
#' @export
boot_pinecone <- function(tree_list, thresh, rthreshold, quiet=TRUE){

  # Some checks
  if(!class(tree_list) %in% c("list", "multiPhylo")) stop("tree_list must be a list or 'multiPhylo' of phylo objects!")
  if(!all(unlist(lapply(tree_list, function(tree) class(tree)=="phylo")))) stop("tree_list must be a list or 'multiPhylo' of phylo objects!")
  tip.lables <- tree_list[[1]]$tip.label
  if(!all(unlist(lapply(tree_list, function(tree) {
    return(all(tree$tip.label %in% tip.lables) && (length(tree$tip.label)==length(tip.lables)))}
  )))) stop("all trees must contain the same taxa!")
  if(!(is.numeric(thresh) && thresh>0)) stop("thresh must be a positive integer!")
  if(!(is.numeric(rthreshold) && rthreshold>0)) stop("rthreshold must be a positive integer!")

  # run pinecone on each phylogeny
  pine_results <- lapply(tree_list, function(tree) {
    return(pinecone(tree, thresh, rthreshold, quiet)$table)
  })

  # summarise into matrix
  occ_matrix <- matrix(0, ncol=length(tip.lables), nrow = length(tip.lables),
                       dimnames = list(tip.lables, tip.lables))
  for (clusters in pine_results){
    for(clust in unique(clusters[,2])){
      if (sum(clusters[,2]==clust)>1){
        pw <- t(combinat::combn(clusters[clusters[,2]==clust,1], 2))
        occ_matrix[pw] <- occ_matrix[pw] + 1
        occ_matrix[pw[,c(2,1),drop=FALSE]] <- occ_matrix[pw[,c(2,1),drop=FALSE]] + 1
      }
    }
  }
  diag(occ_matrix) <- length(tree_list)

  return(occ_matrix)
}
