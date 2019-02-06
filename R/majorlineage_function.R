#' Major lineage
#'
#' Function to identify major lineages composed of identified lineages
#' Singletons which may not have a lineage can also be inlcuded in
#' major lineage
#'
#'@param slnum Number of Sub-lineage identified.
#'@param tree The rooted phylogenetic tree.
#'@param assign Vector of lineage numbers corresponding to tips on a tree.
#'@param ntips Number of tips on the phylogenetic tree.
#'@param rthreshold Number of internal nodes to the root each lineage must have
#' intersecting to be declared as a major sub-lineage.
#'@param treeNnodes Number of internal nodes on phylogenetic tree.
#'
#'majorlineage()

#===========================================================================#
#                                                                           #
#                   Identification of Major Sub-Lineages                    #
#                                                                           #
#===========================================================================#

majorlineage <- function(slnum, tree, assign, rthreshold){

  ntips <- Ntip(tree)
  treeNnodes <- tree$Nnode
  sl_intersect_list <- list()

  if(slnum<=1){
    # we cant find any major sublineages -> output variables
    output <- list(
      maj_sublineage_assign = rep(NA, ntips),
      majorsublist = NA
    )
    return(output)
  }

  for (i in 1:slnum){
    # For each SG retrieve relevant indexes
    slnum_elements <- which(assign == i)
    slnum_elements_tips <- slnum_elements[which(slnum_elements < ntips + 1)]
    iteration <- 0
    sl_element_ancestors_list <- list()

    # Retrieve the internal nodes of each member of said SG
    for (x in slnum_elements_tips){

      # Generate a list of each SG element
      sl_element_ancestors <- Ancestors(tree, x, type = c("all"))
      iteration <- iteration + 1

      # Store the list of ancestors
      sl_element_ancestors_list[iteration] <- list(sl_element_ancestors)
    }

    # Collapse the list into intersects and store for that slnum
    sub_intersect <- Reduce(intersect, sl_element_ancestors_list)

    # Store the list of intersecting ancestors
    sl_intersect_list[i] <- list(sub_intersect)
  }

  # Combinations in columns
  combinations <- combn(length(sl_intersect_list), 2)

  # Take the ancestor nodes of each lineage and setup a paired-wise comparison  - list of lists of the combinations
  ll <- combn(sl_intersect_list, 2, simplify = FALSE )

  # Pair-wise Comparison Function! Intersect the list elements - find the intersect &  the length
  out <- lapply( ll, function(x) length( intersect( x[[1]], x[[2]] ) ) )

  # Declare an empty variable to store wanted comparisons
  major_sublineage_intersect <- list()

  # which lists in the Pair-wise comparison list of lists have equal or more than internal node relatibility threshold
  major_sublineage_comb <- which(out >= rthreshold)

  if(length(major_sublineage_comb)<=0){
    # we cant find any major sublineages -> output variables
    output <- list(
      maj_sublineage_assign = rep(NA, ntips),
      majorsublist = NA
    )
    return(output)
  }

  # Extracting the subgroup number for the above
  if(length(major_sublineage_comb)>0){
    for (i in 1:length(major_sublineage_comb)){

      combin_element <- combinations[, major_sublineage_comb[i]]
      major_sublineage_intersect[[i]] <- combin_element
    }
  }

  # Combination Function! need to combine lists with overlapping slnum
  # Goal: If there are overlapping numbers in differents lists combine them
  for (i in seq_along(major_sublineage_intersect)
       [-length(major_sublineage_intersect)]){
    if (length(intersect(major_sublineage_intersect[[i]],
                         major_sublineage_intersect[[i + 1]])) > 0) {

      major_sublineage_intersect[[i + 1]] <- sort.int(unique(c(major_sublineage_intersect[[i]],
                                                             major_sublineage_intersect[[i + 1 ]])))
      major_sublineage_intersect[[i]] <- as.list(NULL)

    }
  }

  # Take only the rows with something in it - The major subgroup list is then created
  major_sublineage <- Filter(function(x) length(x) > 0, major_sublineage_intersect)

  # Retrieving the SubGroup numbers for each Major SubGroup,
  # Retrieve the SubGroup's tips and storing them in a vector according to the tree
  maj_sublineage_assign <- rep(0, ntips + treeNnodes)

  for (x in 1:length(major_sublineage)){

    maj_sublineage_mems <- NULL

    # Within the Major lineage, call the tips of each lineage and store them in "maj_sublineage_mems"
    for (i in major_sublineage[[x]]){
      maj_sublineage_mems <- c(maj_sublineage_mems, which(assign == i))
    }

    # Similar to the "assign" vector, for each tip noted with the major lineage in the positions of the tips
    maj_sublineage_assign[maj_sublineage_mems] <- x

  }

  # Recording the slnums and their respective major lineage for future reference
  slnum_major_sublineage_list <- matrix(ncol = 2, nrow = length(unlist(major_sublineage)))
  slnum_major_sublineage_list_interation <- 0

  # Each list is a major subgroup i.e. [[1]] is the 1st major subgroup
  for (i in 1:length(major_sublineage)){

    # Extracting each slnum
    for (j in major_sublineage[[i]]){

      slnum_major_sublineage_list_interation <- slnum_major_sublineage_list_interation + 1

      # Appending each row with slnum and respective major lineage
      slnum_major_sublineage_list[slnum_major_sublineage_list_interation, ] <- c(j, i)

    }

  }

  #===========================================================================#
  #                                                                           #
  #       Identification of Major Lineages: Singleton Analysis                #
  #                                                                           #
  #===========================================================================#


  # Check all singletons if they are related to isolates in a lineage;
  # If that lineage is apart of the Major lineage, assign singleton to Major lineage

  # Step 1: Retrieve ancestor nodes of each singleton and store them
  singleton_ancestor_list <- NULL

  # Singletons are those without an assign slnum
  singleton_elements <- which(assign[1:ntips] == 0)

  # Declare how may singletons
  nsingleton <- length(singleton_elements)

  # Combinations in columns
  combinations <- combn(length(sl_intersect_list), 2)

  # Retrieve ancestors for  said singleton & store into a list
  if (nsingleton>0){
    for (i in 1:nsingleton){
      said_singleton <- singleton_elements[i]

      # Generate a list of each slnum taking the 1st element
      singleton_ancestors <- Ancestors(tree, said_singleton, type = c("all"))

      # List in order of singleton number
      singleton_ancestor_list[i] <- list(singleton_ancestors)
    }
  }

  # Step 2: Create a basis to compare the ancestor nodes of each singleton with those of each lineage

  result.df <- expand.grid(sl_intersect_list, singleton_ancestor_list)

  # For each singleton, retrieve the "block" of said singleton + lineage number combinations
  if (nsingleton>0){
    for (i in 1:nsingleton){

      # Define the block
      singleton_slnum_length <- rep(0, slnum)
      single_block_start <- (slnum * i) - slnum + 1
      single_block_end <- i * slnum
      single_block <- result.df[single_block_start:single_block_end, ]

      # x is the slnum number
      for (x in 1:nrow(single_block)){

        # Determine the intersect for a said row
        int <- Reduce(intersect, lapply(single_block, "[[", x))

        # Storing how many internal nodes of a lineage intersect with said singleton over or equal to set Relibility Threshold
        if (length(int) >= rthreshold){

          singleton_slnum_length[[x]] <- length(int)

        }

        # Which lineages overlap in internal nodes with said singleton no equal or more to set threshold
        slnum_for_singleton_int <- which(singleton_slnum_length >= rthreshold)

        # Singleton needs to have intersected with a lineage && the lineage needs to be declared under a Major lineage
        if (length(slnum_for_singleton_int) > 0 && any(slnum_major_sublineage_list[, 1] %in% slnum_for_singleton_int) ){

          # Retrieve corresponding Major lineage
          majsl_of_singleton <- slnum_major_sublineage_list[ which(slnum_major_sublineage_list[, 1] %in% slnum_for_singleton_int), 2]

          # If a singleton has a unique major lineage, store it in the major lineage assigning vector
          # results are a vector calling a Major lineage multiple times i.e. singleton can have same ancestral nodes with multiple slnums

          # check if major lineage is unique, should have only one calling of a major lineage number.
          if (length(unique(majsl_of_singleton)) == 1){

            # assign the major lineage to the singletons
            majlineage <- unique(majsl_of_singleton)

            # convert singleton number into the actual singleton tip number
            singleton_for_assign <- singleton_elements[i]

            # assign the major lineage identified to the singleton
            maj_sublineage_assign[singleton_for_assign] <- majlineage

          } else {
            stop(paste("Error: Two Major Subgroups identified for a singleton ", unique(majsl_of_singleton)), sep = "\n")
          }

        } # End of If Statement: Major lineage of a Singleton Identification

      }
    }

  }


  #output variables
  output <- list(

    maj_sublineage_assign = maj_sublineage_assign,

    majorsublist = slnum_major_sublineage_list

  )

  return(output)
}
