print("testthat.R running")
library(utils)
library(testthat)
library(maria)
if (!require(geiger))  {
  install.packages("geiger", repos = "http://www.stats.bris.ac.uk/R/")
  library("geiger") # ape and gieger to read the tree
}

if (!require(ape)) {
  install.packages("ape", repos = "http://www.stats.bris.ac.uk/R/")
  library("ape")
}

if (!require(igraph)) {
  install.packages("igraph", repos = "http://www.stats.bris.ac.uk/R/")
  library("igraph")
} # igraph for the search

if (!require(BMhyd)) {
  install.packages("BMhyd", repos = "http://www.stats.bris.ac.uk/R/")
  library("BMhyd")
} # loaded to use 'GetAncestor' function

if (!require(phangorn)) {
  install.packages("phangorn", repos = "http://www.stats.bris.ac.uk/R/")
  library("phangorn")
}

if (!require(phytools)) {
  install.packages("phytools", repos = "http://www.stats.bris.ac.uk/R/")
  library("phytools")
}
load_tree(system.file("extdata", "itol_rooted_SNPsites_map_1036491_341_G431_rerun_131117_acctran_steps_211117.tre", package = "testdat"))
sample_id_file <- system.file("extdata", "sample_id.csv", package = "testdat")

thresh <- maria(tree,2,2,sample_id_file)

# test_check("maria")
# test_check("get_path_distance_function")
# test_check("load_tree_function")

test_that("Input for thresh is recognised", {
  expect_equal(thresh == 2)
})
