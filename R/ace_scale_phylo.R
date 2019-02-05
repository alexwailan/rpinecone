#' ace_scale_phylo
#'
#' A wrapper around phangorn functions to scale the branch lengths of a phylogeny by the number of loci seperating ancestral states.
#'
#'
#' @param tree the phylogeny to be scaled
#' @param data the sequence data associated with the same phylogeny
#'
#' @return A phylo object with branch lengths scaled by SNP distance.
#'
#' @examples
#'
#'
#' @export
ace_scale_phylo <- function(tree, data){

  # check inputs
  if(class(tree)!="phylo") stop("tree is not a phylo object!")
  if(class(data)!="DNAbin") stop("data is not a DNAbin object!")

  # ancestral state reconstruction
  phydata <- phangorn::phyDat(data)
  anc_states <- phangorn::ancestral.pars(tree, phydata, type = "ACCTRAN")

  # adjust edge lengths
  tree$edge.length <- unlist(lapply(1:nrow(tree$edge), function(i){
    sum(anc_states[[tree$edge[i,1]]]!=anc_states[[tree$edge[i,2]]])/2
  }))

  return(tree)
}
