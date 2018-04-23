# pinecone

This package defines sub-lineages within a bacterial clonal expansion via a phylogenetic tree.

Inputs for the package 
	- A rooted accelerated transformation (ACCTRANS) tree in newick format
	- SNP threshold - used for sub-grouping
	- Relatibility threshold - the number of internal vertices each sub-lineage must have between each other to form a Major Sub-group


# Output of pinecone
Output is a list of lists:

	$SLNo: Number of Subgroups Identified

	$majorSLno: Number of Major Subgroups

	$singletons: Number of Singletons remaining

	$table: Table of the subgrouping analysis stating the Sub-group and Major Sub-group identified for each isolate.

	$ntips: The number of tips in the tree

# Output for iTOL
Using the output from pinecone, there are two functions which can parse the output of the analysis to display on a tree on iTOL. Below is an example and description.

	output <- pinecone(tree,5,3)

	itol_labels_template(output) - Exports a data file for replacing Tip Labels with Sub-Group number (LABELS)
	
	itol_sublineage_output(output) - Exports a data file for displaying the Sub-Groups (DATASET_COLOURSTRIP)

	itol_major_SL_output(output) -  Exports a data file for displaying the Major Sub-Groups (DATASET_COLOURSTRIP)

	