#' Major Subs
#'
#' Function to identify major sub-groups composed of identified sub-groups
#' Singletons which may not have a sub-group can also be inlcuded in
#' major sub-group
#'
#'@param sgnum Number of Sub-groups identified.
#'@param tree The phylogenetic tree.
#'@param assign Vector of Sub-group numbers corresponding to tips on a tree.
#'@param ntips Number of tips on the phylogenetic tree.
#'@param rthreshold Number of internal nodes to the root each sub-group must have
#' intersecting to be declared as a major sub-group.
#'@param treeNnodes Number of internal nodes on phylogenetic tree.
#'
#' majorsub()

#===========================================================================#
#                                                                           #
#                   Identification of Major Sub-Grouping                    #
#                                                                           #
#===========================================================================#

majorsub <- function(sgnum, tree, assign, ntips, rthreshold, treeNnodes){

  #treeNnodes is ode - need to fix
  sg_intersect_list <- list()

  for (i in 1:sgnum){
    #for each SG retrieve relevant indexes
    sgnum_elements <- which(assign == i)
    sgnum_elements_tips <- sgnum_elements[which(sgnum_elements < ntips + 1)]
    iteration <- 0
    sg_element_ancestors_list <- list()

    #Retrieve the internal nodes of each member of said SG
    for (x in sgnum_elements_tips){

      #generate a list of each SG element
      sg_element_ancestors <- Ancestors(tree, x, type = c("all"))
      iteration <- iteration + 1

      #Store the list of ancestors
      sg_element_ancestors_list[iteration] <- list(sg_element_ancestors)
    }

    #Collapse the list into intersects and store for that sgnum
    sub_intersect <- Reduce(intersect, sg_element_ancestors_list)

    #Store the list of intersecting ancestors
    sg_intersect_list[i] <- list(sub_intersect)
  }

  # Combinations in columns
  combinations <- combn(length(sg_intersect_list), 2)

  # Take the ancestor nodes of each Sub-Group and setup a paired-wise comparison  - list of lists of the combinations
  ll <- combn(sg_intersect_list, 2, simplify = FALSE )

  # Pair-wise Comparison Function! Intersect the list elements - find the intersect &  the length
  out <- lapply( ll, function(x) length( intersect( x[[1]], x[[2]] ) ) )

  # Declare an empty variable to store wanted comparisons
  major_subgroup_intersect <- list()

  # which lists in the Pair-wise comparison list of lists have equal or more than internal node relatibility threshold
  major_subgroup_comb <- which(out >= rthreshold)

  # Extracting the subgroup number for the above
  for (i in 1:length(major_subgroup_comb)){

    combin_element <- combinations[, major_subgroup_comb[i]]
    major_subgroup_intersect[[i]] <- combin_element

  }

  # Combination Function! need to combine lists with overlapping sgnum;
  # Goal: If there are overlapping numbers in differents lists combine them
  for (i in seq_along(major_subgroup_intersect)
       [-length(major_subgroup_intersect)]){
    if (length(intersect(major_subgroup_intersect[[i]],
                         major_subgroup_intersect[[i + 1]])) > 0) {

      major_subgroup_intersect[[i + 1]] <- sort.int(unique(c(major_subgroup_intersect[[i]],
                                                             major_subgroup_intersect[[i + 1 ]])))
      major_subgroup_intersect[[i]] <- as.list(NULL)

    }
  }

  # Take only the rows with something in it - The major subgroup list is then created
  major_subgroup <- Filter(function(x) length(x) > 0, major_subgroup_intersect)

  # Retrieving the SubGroup numbers for each Major SubGroup,
  # Retrieve the SubGroup's tips and storing them in a vector according to the tree
  maj_subgroup_assign <- rep(0, ntips + treeNnodes)

  for (x in 1:length(major_subgroup)){

    maj_subgroup_mems <- NULL

    # Within the Major Sub-Group, call the tips of each Sub-Group and store them in "maj_subgroup_mems"
    for (i in major_subgroup[[x]]){
      maj_subgroup_mems <- c(maj_subgroup_mems, which(assign == i))
    }

    # Similar to the "assign" vector, for each tip noted with the major cluster in the positions of the tips
    maj_subgroup_assign[maj_subgroup_mems] <- x

  }

  # Recording the SGnums and their respective major Sub-Group for future reference
  sgnum_major_subgroup_list <- matrix(ncol = 2, nrow = length(unlist(major_subgroup)))
  sgnum_major_subgroup_list_interation <- 0

  # each list is a major subgroup i.e. [[1]] is the 1st major subgroup
  for (i in 1:length(major_subgroup)){

    # extracting each SGnum
    for (j in major_subgroup[[i]]){

      sgnum_major_subgroup_list_interation <- sgnum_major_subgroup_list_interation + 1

      # appending each row with sgnum and respective major cluster
      sgnum_major_subgroup_list[sgnum_major_subgroup_list_interation, ] <- c(j, i)

    }

  }

  #===========================================================================#
  #                                                                           #
  #       Identification of Major Sub-Grouping: Singleton Analysis            #
  #                                                                           #
  #===========================================================================#


  # Check all singletons if they are related to isolates in a Sub-Group;
  #  if that cluster is apart of the Major Sub-Group, assign singleton to Major Sub-Group

  # Step 1: Retrieve ancestor nodes of each singleton and store them
  singleton_ancestor_list <- NULL

  # Singletons are those without an assign sgnum
  singleton_elements <- which(assign[1:ntips] == 0)

  # Declare how may singletons
  nsingleton <- length(singleton_elements)

  # combinations in columns
  combinations <- combn(length(sg_intersect_list), 2)

  # Retrieve ancestors for  said singleton & store into a list
  for (i in 1:nsingleton){
    said_singleton <- singleton_elements[i]

    # generate a list of each sgnum taking the 1st element
    singleton_ancestors <- Ancestors(tree, said_singleton, type = c("all"))

    # list in order of singleton number
    singleton_ancestor_list[i] <- list(singleton_ancestors)

  }

  # Step 2: Create a basis to compare the ancestor nodes of each singleton with those of each Sub-group

  result.df <- expand.grid(sg_intersect_list, singleton_ancestor_list)

  # for each singleton retrieve the "block" of said singleton + sgnum combinations
  for (i in 1:nsingleton){

    # Define the block
    singleton_sgnum_length <- rep(0, sgnum)
    single_block_start <- (sgnum * i) - sgnum + 1
    single_block_end <- i * sgnum
    single_block <- result.df[single_block_start:single_block_end, ]

    # x is the sgnum number
    for (x in 1:nrow(single_block)){

      # Determine the intersect for a said row
      int <- Reduce(intersect, lapply(single_block, "[[", x))

      # Storing how many internal nodes of a Sub-Group intersect with said singleton over or equal to set Relibility Threshold
      if (length(int) >= rthreshold){

        singleton_sgnum_length[[x]] <- length(int)

      }

      # which Sub-groups overlap in internal nodes with said singleton no equal or more to set threshold
      sgnum_for_singleton_int <- which(singleton_sgnum_length >= rthreshold)

      # singleton needs to have intersected with a Sub-Group && the Sub-Group needs to be declared under a Major Sub-Group
      if (length(sgnum_for_singleton_int) > 0 && any(sgnum_major_subgroup_list[, 1] %in% sgnum_for_singleton_int) ){

        # Retrieve corresponding Major Sub-group
        majsb_of_singleton <- sgnum_major_subgroup_list[ which(sgnum_major_subgroup_list[, 1] %in% sgnum_for_singleton_int), 2]

        # If a singleton has a unique major cluster, store it in the major cluster assigning vector
        # results are a vector calling a Major cluster multiple times i.e. singleton can have same ancestral nodes with multiple sgnums

        # check if major cluster is unique, should have only one calling of a major cluster number.
        if (length(unique(majsb_of_singleton)) == 1){

          # assign the major cluster to the singletons
          majclust <- unique(majsb_of_singleton)

          # convert singleton number into the actual singleton tip number
          singleton_for_assign <- singleton_elements[i]

          # assign the major cluster identified to the singleton
          maj_subgroup_assign[singleton_for_assign] <- majclust

        } else {
          cat(paste("Error: Two Major Subgroups identified for a singleton ", unique(majsb_of_singleton)), sep = "\n")
        }

      } # End of If Statement: Major Cluster of a Singleton Identification

    }

  }
  output <- list(

    maj_subgroup_assign = maj_subgroup_assign,

    majorsublist = sgnum_major_subgroup_list

  )

  return(output)
}
