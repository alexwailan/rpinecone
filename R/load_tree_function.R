#' Load a Newick tree
#'
#' This function loads a phylogenetic tree. The tree must be rooted, newick
#' format and an ACCTRANS tree.
#' @import phytools
#' @param infile Path to the input file
#' @return A matrix of the infile

load_tree <- function(infile){
  tree <- read.newick(infile)
  return(tree)
}
