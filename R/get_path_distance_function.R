#' Get the distance of a path
#'
#' This function returns the path distance between node and tip on a rooted tree, in a root to tip direction.
#' @param graph_tree A rooted tree of newick format, converted into graph form by igraph.
#' @param root A node in the graph you want to begin the path.
#' @param node target node.
#' @param input_tree a rooted tree that has been parsed into this script - class "phylo".
#' @import igraph
#' @return The SNP distance between given root and node
#' get_path_distance()

get_path_distance <- function(graph_tree, root, node, input_tree){
  path <- shortest_paths(graph_tree,
                         root,
                         node,
                         mode = "out",
                         weights = NULL,
                         output = "both") #Determine path between nodes.

  #retrieve edgelist of tree in igraph form (different from the original tree)
  igraph_edgelist <- as_edgelist(graph_tree)

  #define a matrix to store edges
  edge_matrix <-  matrix(ncol = 2, nrow = length(path$epath[[1]]))

  # retrieve edges of path and store them
  for (i in 1:length(path$epath[[1]])){

    #stating the edge; retrieved as "row" number in edgelist vector
    edge <- as.integer(path$epath[[1]][i])

    #retrieving the edge from list
    edge_matrix[i, ] <- igraph_edgelist[edge, ]
  }
  #comparing two rows; if both elements match i.e. all() will call TRUE
  step1 <- Vectorize(function(x, y) {
    all(input_tree$edge[x, ] == edge_matrix[y, ])
  })

  #Generate a matrix of a pairwise comparison between two matrices and their edges
  #Rows correspond to tree$edges vs cols correspond to edges in path
  result <- outer(1:nrow(input_tree$edge),
                  1:nrow(edge_matrix),
                  FUN = step1)
  #identify which row each edge is found in input_tree$edge
  edge_row <- which(apply(result, 1, any))

  path_length <- 0

  #for each edge row, return its corresponding edge length and add to the path length
  for (i in edge_row){
    edge_length <- as.integer(input_tree$edge.length[i])
    path_length <- path_length + edge_length
  }

  return(path_length)

} #End of Function
