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

  # Resolving dichotomies by deleting all branches smaller than 0.5 i.e. zero SNPs and collapses the corresponding dichotomies into a multichotomy.
  tree <- di2multi(tree, 0.5)
  tree_name <- paste("NewickTreefile_", Date, ".tre")
  write.tree(tree, file = tree_name) # Output the new tree.

  # Display Number of Taxa in Tree.
  cat("\t", length(tree$tip.label), " taxa found from phylogenetic tree\n", sep = "")

  # Renaming the node labels; at a sequence from 1 to (number of nodes in tree), at increments of 1.
  tree$node.label <- paste("Node", seq(1, tree$Nnode, 1), sep = "_")

  #===========================================================================#
  #                                                                           #
  #                  Define variables for Sub-grouping                        #
  #                                                                           #
  #===========================================================================#

  # Number of tips
  ntips <- Ntip(tree)

  # setting a Sub-Group number - will tick up as more clusters are identified
  sgnum <- 0

  # cluster assignment - rep replicates the values found in x. rep(x, ...);
  # generating a vector the size of ntips+tree$Nnode corresponding to the number of positions within the tree data
  assign <- rep(0, ntips + tree$Nnode)

  # convert tree in igraph form - Takes a 2 column matrix edge list, where each row is an existing edge
  igraph.tree <- graph.edgelist(tree$edge)

  # Depth First Search traverses a graph, begins at a "root" vertex and tries to go quickly as far from as possible
  # Order - return the DFS ordering of the vertices; dist - return the distance from the root of the search tree
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

          # Some internal nodes have downstream tips with a distance
          # over declared SNP threshold,
          # however have tips in zero distance to node
          # Assign these tips as  a Sub-Group

          sgnum <- sgnum + 1

          assign[tips_0_dist_node] <- sgnum

          assign[node] <- sgnum

          cat(paste("Internal sub-group detected - Sub-group number assigned: ", sgnum), sep = "\n")

        }

      } #End of If statement

  } # End of for loop

  #===========================================================================#
  #                                                                           #
  #                   Identification of Major Sub-Grouping                    #
  #                                                                           #
  #===========================================================================#

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

  # Combination Function! need to combine lists with overlapping sgnum; Goal: If there are overlapping numbers in differents lists combine them
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

  # Number of Major SubGroup
  numbmajorsubgroup <- length(major_subgroup)

  # Retrieving the SubGroup numbers for each Major SubGroup, retrieve the SubGroup's tips and storing them in a vector according to the tree
  maj_subgroup_mems <- NULL
  maj_subgroup_assign <- rep(0, ntips + tree$Nnode)

  for (x in 1:length(major_subgroup)){

    maj_subgroup_mems <- NULL

    # Within the Major Sub-Group, call the tips of each Sub-Group and store them in "maj_subgroup_mems"
    for (i in major_subgroup[[x]]){
     maj_subgroup_mems <- c(maj_subgroup_mems, which(assign == i))
     }

    # Similar to the "assign" vector, for each tip noted with the major cluster in the positions of the tips
    maj_subgroup_assign[maj_subgroup_mems] <- x

  }

  # recording the SGnums and their respective major Sub-Group for future reference
  sgnum_major_subgroup_list <- matrix(ncol = 2, nrow = length(unlist(major_subgroup)))
  sgnum_major_subgroup_list_interation <- 0

  # each list is a major subgroupu i.e. [[1]] is the 1st major subgroup
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

  # Number of singletons after sub-grouping
  remaining_singletons <- which(assign[1:ntips] == 0)

  majorsubgrouptable <- mat.or.vec(length(tree$tip.label), 3)

  # Collating the tip labels and their respective Sub-Group and Major Sub-Group
  colnames(majorsubgrouptable) <- c("Taxa", "Sub-group", "Major.Sub-group")

  majorsubgrouptable[, 1] <- tree$tip.label

  majorsubgrouptable[, 2] <- assign[1:ntips]

  majorsubgrouptable[, 3] <- maj_subgroup_assign[1:ntips]

  majorzero_elems <- which(majorsubgrouptable[, 2] == 0)

  # Labeling Singletons
  majorsubgrouptable[majorzero_elems, 2] <- "singleton_"

  singletonsii <- which(majorsubgrouptable[, 2] == "singleton_")

  majorsubgrouptable[singletonsii, 2] <- paste("singleton_", seq(1, length(singletonsii), 1), sep = "")

   # collating variables from above
   output <- list(

     # 1: Number of Major Subgroups
     majorSBno = numbmajorsubgroup,

     # 2: Singleton nodes
     singletons = remaining_singletons,

     # 3: Itol output table
     itolOutput = majorsubgrouptable,

     # 4: stating the number of tips in the tree
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
     numbmajorsubgroupmems <- length(which(majorsubgrouptable[, 3] == i))
     numbmajorsubgroupsgnum <- length(which(sgnum_major_subgroup_list[, 2] == i))
     cat(paste("Major Sub-Group ", i, "is composed of ", numbmajorsubgroupsgnum, " Sub-Groups & ", numbmajorsubgroupmems, " isolates."), sep = "\n")
   }

   return(output)

}
