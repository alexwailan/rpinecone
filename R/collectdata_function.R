#'collectiondata
#'
#'Collating data for the output.
#'
#'@param tree The rooted phylogenetic tree.
#'@param assign Vector of Sub-group numbers corresponding to tips on a tree.
#'@param maj_sublineage The output of the Major Sub-Lineage function
#'
#' collectiondata()


collectdata <- function(tree,assign,maj_sublineage){

  table <- mat.or.vec(length(tree$tip.label), 3)
  ntips <- Ntip(tree)

  # Collating the tip labels and their respective Sub-Group and Major Sub-Group
  colnames(table) <- c("Taxa", "Sub-lineage", "Major.Sub-lineage")

  table[, 1] <- tree$tip.label

  table[, 2] <- assign[1:ntips]

  table[, 3] <- maj_sublineage$maj_sublineage_assign[1:ntips]

  elems <- which(table[, 2] == 0)

  # Labeling Singletons
  table[elems, 2] <- "singleton_"

  singletonsii <- which(table[, 2] == "singleton_")

  if(length(singletonsii)>0){
    table[singletonsii, 2] <- paste("singleton_", seq(1, length(singletonsii), 1), sep = "")
  }

  return(table)
}
