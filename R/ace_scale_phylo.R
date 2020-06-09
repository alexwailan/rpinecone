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

  # convert all non acgt's to unknowns
  for (i in 1:nrow(data)){
    data[i,!(data[i,] %in% as.DNAbin(c('a','c','g','t')))] <- as.DNAbin('-')
  }

  # ancestral state reconstruction
  phydata <- phangorn::phyDat(data)
  anc_states <- phangorn::ancestral.pars(tree, phydata, type = "ACCTRAN")

  # scale edge lengths to snp distance
  tree$edge.length <- unlist(lapply(1:nrow(tree$edge), function(i){
    sum((anc_states[[tree$edge[i,1]]]!=anc_states[[tree$edge[i,2]]]) &
          (anc_states[[tree$edge[i,1]]]!=0.25) &
          (anc_states[[tree$edge[i,2]]]!=0.25))/2
  }))

  return(tree)
}
