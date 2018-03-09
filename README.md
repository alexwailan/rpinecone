# maria

This package sub-groups bacterial isolates of a clonal expansion via a phylogenetic tree.

Inputs for the package 
	- A rooted accelerated transformation (ACCTRANS) tree in newick format
	- SNP threshold - used for sub-grouping
	- Relatibility threshold - the number of internal nodes to the root each sub-group must have between each other to be declared as a Major Sub-group


# Output of pinecone
Output is a list of lists:

	$SGNo: Number of Subgroups Identified

	$majorSBno: Number of Major Subgroups

	$singletons: Number of Singletons remaining

	$output: Table of the subgrouping analysis stating the Sub-group and Major Sub-group identified for each isolate.

	$ntips: The number of tips in the tree

# Output for iTOL
Using the output from pinecone, there are two functions which can parse the output of the analysis to display on a tree on iTOL. Below is an example and description.

	output <- pinecone(tree,5,3)

	itol_labels_template(output) - Exports a data file for replacing Tip Labels with Sub-Group number (LABELS)

	itol_major_SB_binary_template(output) -  Exports a data file for displaying the Major Sub-Groups (DATASET_COLOURSTRIP)

	itol_subGroup_binary_template(output) - Exports a data file for displaying the Sub-Groups (DATASET_COLOURSTRIP)