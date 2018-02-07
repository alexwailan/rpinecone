library(utils)
library(testthat)
library(maria)

library(geiger) # ape and gieger to read the tree

library(ape)
library(igraph)
library(BMhyd)
# loaded to use 'GetAncestor' function
library(phangorn)
library(phytools)

expect_true (file.exists (file.path ("inst",
                                     "extdata",
                                     "itol_rooted_SNPsites_map_1036491_341_G431_rerun_131117_acctran_steps_211117.tre")))

tree <- file.path("inst",
                            "extdata",
                            "itol_rooted_SNPsites_map_1036491_341_G431_rerun_131117_acctran_steps_211117.tre")



sample_id_file <- file.path("inst",
                            "extdata", "sample_id.csv", package = "testdat")

maria(tree,2,2)

print(thresh)
# test_check("maria")
# test_check("get_path_distance_function")
# test_check("load_tree_function")


