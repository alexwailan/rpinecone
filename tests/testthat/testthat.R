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
library(RColorBrewer)

print (file.path())
expect_true (file.exists
             (file.path
               ("..",
               "..",
               "inst",
               "extdata",
               "itol_rooted_SNPsites_map_1036491_341_G431_rerun_131117_acctran_steps_211117.tre")))

tree <- file.path(#"..",
                  #"..",
                  "inst",
                  "extdata",
                  "itol_rooted_SNPsites_map_1036491_341_G431_rerun_131117_acctran_steps_211117.tre")



sample_id_file <- file.path("inst",
                            "extdata",
                            "sample_id.csv",
                            package = "testdat")

output <- maria(tree,2,3)

itol_labels_template(output$ntips,output$itolOutput)

itol_major_SB_binary_template(output$majorSBno,output$ntips,output$itolOutput)

# test_check("maria")
# test_check("get_path_distance_function")
# test_check("load_tree_function")


