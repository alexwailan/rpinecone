#' maria
#' @param tree Newick formatted ACCTRANS tree.
#' @param thresh SNP threshold for sub-grouping isolates.
#' @param rthreshold The threhold number of ancestral nodes to the root
#' each sub-group must have if they are to be declared as major sub-groups
#' @param sample_id_file CSV file, with a list of sample ids in the first column
#' and sample ids in the second column
#' @import geiger ape igraph BMhyd phangorn phytools

maria <- function(tree,thresh,rthreshold){

  #Retrieving date for output file naming purposes.
  Date <- Sys.Date()

  tree <- read.newick(tree)


  #Resolving dichotomies by deleting all branches smaller than 0.5 i.e. zero SNPs and collapses the corresponding dichotomies into a multichotomy.
  tree <- di2multi(tree, 0.5)
  tree_name <- paste("NewickTreefile_", Date, ".tre")
  write.tree(tree, file = tree_name) #Output the new tree.

  #Display Number of Taxa in Tree.
  cat("\t", length(tree$tip.label), " taxa found from phylogenetic tree\n", sep = "")

  #Renaming the node labels; at a sequence from 1 to (number of nodes in tree), at increments of 1.
  tree$node.label <- paste("Node", seq(1, tree$Nnode, 1), sep = "_")


  #===============================================================================================#
  #                                                                                               #
  #                               Define variables for Sub-grouping                               #
  #                                                                                               #
  #===============================================================================================#

  #Number of tips
  ntips <- Ntip(tree)

  # setting a Sub-Group number - will tick up as more clusters are identified
  sgnum <- 0

  # cluster assignment - rep replicates the values found in x. rep(x, ...); generating a vector the size of ntips+tree$Nnode corresponding to the number of positions within the tree data
  assign <- rep(0, ntips + tree$Nnode)

  # convert tree in igraph form - Takes a 2 column matrix edge list, where each row is an existing edge
  igraph.tree <- graph.edgelist(tree$edge)

  #Depth First Search traverses a graph, begins at a "root" vertex and tries to go quickly as far from as possible
  #Order - return the DFS ordering of the vertices; dist - return the distance from the root of the search tree
  dfs <- graph.dfs(igraph.tree,
                   root = ntips + 1,
                   neimode = "out",
                   order = TRUE,
                   dist = TRUE)

  #===============================================================================================#
  #                                                                                               #
  #                                  Investigating edge distances                                 #
  #                                                                                               #
  #===============================================================================================#

  #Traverse the tree in depth first order - starting at the root
    for (i in 1:length(dfs$order)){
      node <- dfs$order[i]

      # Skip leaves/tips
      #If the number of the "node" is less than ntips+1 it is actually a leaf/tip
      #i.e. phytools stores tips of trees as 1:n (n = no. of tips) then nodes as n+1:n+m (m = no. of nodes)
      if (node < ntips + 1) {
        cat(paste(tree$tip.label[node], " skipped as this is a leaf."), sep = "\n")
         next
      } else {
        #if the number of the "node" is equal or more than ntips+1 it is actually a node; proceed with below
        cat("", paste(tree$node.label[node - ntips], " will be investigated."), sep = "\n")
      }

      #Assign a subgroup number containing the said node & all the tips connected in the 'out' direction from the root
      if (assign[node] <= 0){

        #Said node becomes the root; dfs will search out and notes the order of "nodes"
        #(which includes internal nodes + tips); display only the $order and store it
        subtree <- graph.dfs(igraph.tree, node, neimode = "out", unreachable = FALSE)$order

        #Checks within the igraph.vs called subtree which values do not equal to NA and stores them under subtree i.e. remove all NAs
        subtree <- subtree[!is.na(subtree)]

        #Define maximum distance amongst all tips from said node
        node_tips_max_dist <- node_tips_max_dist(node, subtree, igraph.tree, tree)

        #Define the tips connected to node at zero distance
        tips_0_dist_node <- tips_zero_distance_node(node, subtree, igraph.tree, tree)

        #Logical sub-grouping numbering; composed of  two parameters - threshold && zero distance

        #Logical: tip distances from node below or equal to threshold?
        if (node_tips_max_dist <= thresh){

          #Tick up the clustering number
          sgnum <- sgnum + 1

          #Record/assign the cluster number to the "assign" vector in the positions specified by the "subtree" variable
          assign[subtree] <- sgnum

          #Display that a Sub-Group has been identified and assigned
          cat(paste("Threhold met - Sub-Group Number assigned: ", sgnum), sep = "\n")

        } else
          if (sum(tips_0_dist_node) > 0 &&
              tips_0_dist_node[1] != sum(tips_0_dist_node)){
          #If said internal node has downstream tips over the declared SNP threshold assess tips directly connected to said node
          #Assign Sub-group if there are tips are zero distance to node; Must not be a singleton

          sgnum <- sgnum + 1

          #assign sgnum to tips at node in zero distance
          assign[tips_0_dist_node] <- sgnum

          #assign sgnum to the node
          assign[node] <- sgnum

          #Display that a Sub-Group has been identified and assigned
          cat(paste("Internal sub-group detected - Sub-group number assigned: ", sgnum), sep = "\n")

        }

      } #End of If statement

  } # End of for loop

  #===============================================================================================#
  #                                                                                               #
  #                            Identification of Major Sub-Grouping                               #
  #                                                                                               #
  #===============================================================================================#

  iteration <- 0
  sg_intersect_list <- list()

  for (i in 1:sgnum){
    #which elements are found for a sgnum
    sgnum_elements <- which(assign == i)
    sgnum_elements_tips <- sgnum_elements[which(sgnum_elements < ntips + 1)]
    iteration <- 0
    sg_element_ancestors_list <- list()

    #Retrieve the internal nodes of each member of said subgroup
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

  #combinations in columns
  combinations <- combn(length(sg_intersect_list), 2)

  #Take the ancestor nodes of each Sub-Group and setup a paired-wise comparison  - list of lists of the combinations
  ll <- combn(sg_intersect_list, 2, simplify = FALSE )

  # Pair-wise Comparison Function! Intersect the list elements - find the intersect &  the length
  out <- lapply( ll, function(x) length( intersect( x[[1]], x[[2]] ) ) )

  #which Pair-wise comparisons have more than two internal nodes

  #Declare an empty variable to store wanted comparisons
  major_subgroup_intersect <- list()

  #which lists in the Pair-wise comparison list of lists have equal or more than internal node relatibility threshold
  major_subgroup_comb <- which(out >= rthreshold)

  #Extracting the subgroup number for the above
  for (i in 1:length(major_subgroup_comb)){

   combin_element <- combinations[, major_subgroup_comb[i]]
   major_subgroup_intersect[[i]] <- combin_element

  }

  #Combination Function! need to combine lists with overlapping sgnum; Goal: If there are overlapping numbers in differents lists combine them
  for (i in seq_along(major_subgroup_intersect)
       [-length(major_subgroup_intersect)]) {
   if (length(intersect(major_subgroup_intersect[[i]],
                        major_subgroup_intersect[[i + 1]])) > 0) {

     major_subgroup_intersect[[i + 1]] <- sort.int(unique(c(major_subgroup_intersect[[i]],
                                                          major_subgroup_intersect[[i + 1 ]])))
     major_subgroup_intersect[[i]] <- as.list(NULL)

   }
  }

  #Take only the rows with something in it - The major subgroup list is then created
  major_subgroup <- Filter(function(x) length(x) > 0, major_subgroup_intersect)

  #Number of Major SubGroup
  numbmajorsubgroup <- length(major_subgroup)

  #Retrieving the SubGroup numbers for each Major SubGroup, retrieve the SubGroup's tips and storing them in a vector according to the tree
  maj_subgroup_mems <- NULL
  maj_subgroup_assign <- rep(0, ntips + tree$Nnode)

  for (x in 1:length(major_subgroup)){

    maj_subgroup_mems <- NULL


    #Within the Major Sub-Group, call the tips of each Sub-Group and store them in "maj_subgroup_mems"
    for (i in major_subgroup[[x]]){
     maj_subgroup_mems <- c(maj_subgroup_mems, which(assign == i))
     }

    #Similar to the "assign" vector, for each tip noted with the major cluster in the positions of the tips
    maj_subgroup_assign[maj_subgroup_mems] <- x

  }

  # recording the SGnums and their respective major Sub-Group for future reference
  sgnum_major_subgroup_list <- matrix(ncol = 2, nrow = length(unlist(major_subgroup)))
  sgnum_major_subgroup_list_interation <- 0

  #each list is a major cluster i.e. [[1]] is the 1st major cluster
  for (i in 1:length(major_subgroup)){

    #extracting each SGnum
    for (j in major_subgroup[[i]]){

      sgnum_major_subgroup_list_interation <- sgnum_major_subgroup_list_interation + 1

      #appending each row with sgnum and respective major cluster
      sgnum_major_subgroup_list[sgnum_major_subgroup_list_interation, ] <- c(j, i)

    }

  }

  #===============================================================================================#
  #                                                                                               #
  #                 Identification of Major Sub-Grouping: Singleton Analysis                      #
  #                                                                                               #
  #===============================================================================================#


  # Check all singletons if they are related to isolates in a Sub-Group; if that cluster is apart of the Major Sub-Group, assign singleton to Major Sub-Group

  #Step 1: Retrieve ancestor nodes of each singleton and store them
  singleton_ancestor_list <- NULL

  #Singletons are those without an assign sgnum
  singleton_elements <- which(assign[1:ntips] == 0)

  #Declare how may singletons
  nsingleton <- length(singleton_elements)

  #combinations in columns
  combinations <- combn(length(sg_intersect_list), 2)

  #Retrieve ancestors for  said singleton & store into a list
  for (i in 1:nsingleton){
   said_singleton <- singleton_elements[i]

   #generate a list of each sgnum taking the 1st element
   singleton_ancestors <- Ancestors(tree, said_singleton, type = c("all"))

   #list in order of singleton number
   singleton_ancestor_list[i] <- list(singleton_ancestors)

  }

  #Step 2: Create a basis to compare the ancestor nodes of each singleton with those of each sgnum

  #Pairwise comparison combinations are based per singleton against all sgnums;
  #ancestor nodes of said singleton will be listed in var2, in combination with the number sgnums and their respective ancestor nodes
  #i.e each singleton comparison is in a set of rows the size the number of sgnums identified
  #e.g. first singleton is row 1 to row(no of cums), second singleton is the next set of rows after to the number of sgnums
  result.df <- expand.grid(sg_intersect_list, singleton_ancestor_list)

  #for each singleton retrieve the "block" of said singleton + sgnum combinations
  for (i in 1:nsingleton){

    #Define the block
    singleton_sgnum_length <- rep(0, sgnum)
    single_block_start <- (sgnum * i) - sgnum + 1
    single_block_end <- i * sgnum
    single_block <- result.df[single_block_start:single_block_end, ]

    #x is the sgnum number
    for (x in 1:nrow(single_block)){

      #Determine the intersect for a said row i.e. Singleton vs sgnum(which is row number)
      int <- Reduce(intersect, lapply(single_block, "[[", x))

      #Storing how many internal nodes of a Sub-Group intersect with said singleton over or equal to set Relibility Threshold
      if (length(int) >= rthreshold){

        singleton_sgnum_length[[x]] <- length(int)

      }

      #which sgnum overlap  equal or more in internal nodes with said singleton no
      sgnum_for_singleton_int <- which(singleton_sgnum_length >= rthreshold)

      #singleton needs to have intersected with a sgnum && the sgnum needs to be have a major cluster to call it
      if (length(sgnum_for_singleton_int) > 0 && any(sgnum_major_subgroup_list[, 1] %in% sgnum_for_singleton_int) ){

        #given the list of each sgnum and their corresponding major Sub-Group - retrieve corresponding major Sub-group
        majsb_of_singleton <- sgnum_major_subgroup_list[
          #rows which desired sgnums
          which(sgnum_major_subgroup_list[, 1] %in% sgnum_for_singleton_int),
          2 #2nd column is where the major cluster is stored
                                                        ]

        #If a singleton has a unique major cluster, store it in the major cluster assigning vector
        #results are a vector calling a Major cluster multiple times i.e. singleton can have same ancestral nodes with multiple sgnums

        #check if major cluster is unique, should have only one calling of a major cluster number.
        if (length(unique(majsb_of_singleton)) == 1){

          #assign the major cluster to the singletons in maj_subgroup_assign
          majclust <- unique(majsb_of_singleton)

          #convert singleton number into the actual singleton tip number
          singleton_for_assign <- singleton_elements[i]

          #assign the major cluster identified to the singleton
          maj_subgroup_assign[singleton_for_assign] <- majclust

        }else{
          cat(paste("Error: Two Major Subgroups identified for a singleton ", unique(majsb_of_singleton)), sep = "\n")
        }

      } #End of If Statement: Major Cluster of a Singleton Identification

    }

  }

   # collating variables from above
   ans <- list(

     #Assigned subgroup number for each leaf/tip
     subgroupmems = assign,

     #sizes for each subgroup; including both leaves and internal nodes
     allcsize = table(assign),

     #sizes for each subgroup; leaves/tips only
     leafclustsize = table(assign[1:ntips]),

     #stating the number of tips in the tree
     ntips = ntips,

     #stating the maximum distance that will define a claded
     threshold = thresh,

     majorsubgroupmems = maj_subgroup_assign)

   #subsetting the "membership" to only include tips/leaves
   ans$subgroupmems <- ans$subgroupmems[1:ntips]
   ans$majorsubgroupmems <- ans$majorsubgroupmems[1:ntips]

   #Number of singletons after sub-grouping
   remaining_singletons <- which(assign[1:ntips] == 0)


  #===============================================================================================#
  #                                                                                                 #
  #                                    Saving subgrouping results                                   #
  #                                                                                                 #
  #===============================================================================================#


  # Compiling all meta data, subgrouping and major sub-groups for output.

  if (exists(opt$sampleid)){
    majorsubgrouptable <- mat.or.vec(length(tree$tip.label), 4)

    #Giving column names
    colnames(majorsubgrouptable) <- c("Taxa", "Sample.ID", "Sub-group", "Major.Sub-group")

    #Assign the 1st column of said matrix or vector to the name of the tips from the tree
    majorsubgrouptable[, 1] <- tree$tip.label

    #matches in majorcluster
    filter_lanes <- match(majorsubgrouptable[, 1], Sample.ID[, 1])

    #reordering sample id via tree order
    Sample.ID.order <- Sample.ID[filter_lanes, ]

    #Assign the 2nd column of said matrix or vector to cluster it has been assigned
    majorsubgrouptable[, 2] <- Sample.ID.order[, 2]

    #Assign the 2nd column of said matrix or vector to cluster it has been assigned
    majorsubgrouptable[, 3] <- ans$subgroupmems

    #Assign the 3rd column of said matrix or vector to Major cluster it has been assigned
    majorsubgrouptable[, 4] <- ans$majorsubgroupmems

    #store positions/elements in the vector are equal to zero i.e. were not assigned a cluster number
    majorzero_elems <- which(majorsubgrouptable[, 3] == 0)

    #Said positions in vector are replaced with "singleton"
    majorsubgrouptable[majorzero_elems, 3] <- "singleton_"

    #store positions/elements in the vector have singleton i.e. were not assigned a cluster number
    singletonsii <- which(majorsubgrouptable[, 3] == "singleton_")

    #rename each position stated in "ii" with singleton plus a unique number using seq(1 to length of ii, by increments of 1)
    majorsubgrouptable[singletonsii, 3] <- paste("singleton_", seq(1, length(singletonsii), 1), sep = "")

  } else {
    majorsubgrouptable <- mat.or.vec(length(tree$tip.label), 3)

    #Giving column names
    colnames(majorsubgrouptable) <- c("Taxa", "Sub-group", "Major.Sub-group")

    #Assign the 1st column of said matrix or vector to the name of the tips from the tree
    majorsubgrouptable[, 1] <- tree$tip.label

    #Assign the 2nd column of said matrix or vector to cluster it has been assigned
    majorsubgrouptable[, 2] <- ans$subgroupmems

    #Assign the 3rd column of said matrix or vector to Major cluster it has been assigned
    majorsubgrouptable[, 3] <- ans$majorsubgroupmems

    #store positions/elements in the vector are equal to zero i.e. were not assigned a cluster number
    majorzero_elems <- which(majorsubgrouptable[, 2] == 0)

    #Said positions in vector are replaced with "singleton"
    majorsubgrouptable[majorzero_elems, 2] <- "singleton_"

    #store positions/elements in the vector have singleton i.e. were not assigned a cluster number
    singletonsii <- which(majorsubgrouptable[, 2] == "singleton_")

    #rename each position stated in "ii" with singleton plus a unique number using seq(1 to length of ii, by increments of 1)
    majorsubgrouptable[singletonsii, 2] <- paste("singleton_", seq(1, length(singletonsii), 1), sep = "")
  }


  majsubgroupoutput <- paste(thresh, "_SNPs_SB_threshold_", numbmajorsubgroup, "_MajSubgroups_", Date, ".txt", sep = "")

  write.table(majorsubgrouptable, file = majsubgroupoutput, sep = ",", col.names = T, row.names = F, quote = F)

  #===============================================================================================#
  #                                                                                                 #
  #                              Display Sub-Grouping stats to terminal                             #
  #                                                                                                 #
  #===============================================================================================#

  cat("", paste("Number of Sub-groups identified: ", sgnum),
      paste("Number of Major Subgroups identified: ", length(major_subgroup)),
      paste("Number of Singletons remain: ", length(remaining_singletons)),
      sep = "\n")

  #Display number of samples in major clusters and cluster in each major cluster
  for (i in 1:numbmajorsubgroup){
    numbmajorsubgroupmems <- length(which(majorsubgrouptable[, 3] == i))
    numbmajorsubgroupsgnum <- length(which(sgnum_major_subgroup_list[, 2] == i))
    cat(paste("Major Sub-Group ", i, "is composed of ", numbmajorsubgroupsgnum, " Sub-Groups & ", numbmajorsubgroupmems, " isolates."), sep = "\n")
  }
}
