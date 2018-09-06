
<!-- README.md is generated from README.Rmd. Please edit that file -->
<img src="inst/vignette-supp/rPineconeLogo_full.png" width="100%" />

rPinecone defines sub-lineages within a bacterial clonal expansion via a phylogenetic tree.

Installation
------------

`rPinecone` is currently available on github. It can be installed with `devtools`

``` r
install.packages("devtools")

devtools::install_github("alexwailan/rpinecone")
```

Quick Start
-----------

Run rPinecone

``` r
# devtools::install_github('alexwailan/rpinecone')
library(rPinecone)
library(phytools)

tree.file.name <- system.file("extdata", "pyjar.staph_ST2371_45_re_itol_7079_1_6_root.joint.tre", 
    package = "rPinecone")
tree <- phytools::read.newick(tree.file.name)
results <- pinecone(tree, 2, 3, quiet = TRUE)
```

rPinecone, is an R package designed to define sub-lineages within closely related LV populations. rPinecone uses root-to-tip directional approach to define sub-lineages within a phylogenetic tree according to SNV distance from the ancestral node.

------------------------------------------------------------------------

Libraries
---------

``` r
library(phytools)
library(rPinecone)
```

The package defines sub-lineages within a bacterial clonal expansion via a phylogenetic tree.

Inputs
------

-   A rooted tree in newick format with branch lengths representing SNP distance

    A SNP scaled tree can be produced using the algorithm of Pupko et al. available in python at [pyjar](https://github.com/simonrharris/pyjar#usage)

-   SNP threshold - used for sub-grouping
-   Relatability threshold - the number of internal vertices each sub-lineage must have between each other to form a Major Sub-group

Operation
---------

After loading a newick formatted tree rPinecone can be called as

``` r
tree.file.name <- system.file("extdata", "pyjar.staph_ST2371_45_re_itol_7079_1_6_root.joint.tre", 
    package = "rPinecone")
tree <- phytools::read.newick(tree.file.name)
results <- pinecone(tree, 2, 3, quiet = TRUE)
```

Exporting results
-----------------

The resulting output can the be parsed for plotting with iTOL

To exports a data file for replacing Tip Labels with Sub-Group number (LABELS) run

``` r
itol_labels_template(results)
```

To export a data file for displaying the Sub-Groups (DATASET\_COLOURSTRIP) run

``` r
itol_sublineage_output(results)
```

To export a data file for displaying the Major Sub-Groups (DATASET\_COLOURSTRIP) run

``` r
itol_major_SL_output(results)
```

The phylogentic tree with dichotomies resolved into multichotomies can be saved by running

``` r
write.tree(results$tree, file = "rpinecone.tree")
```

References
----------

Letunic, Ivica, and Peer Bork. 2007. “Interactive Tree of Life (iTOL): An Online Tool for Phylogenetic Tree Display and Annotation.” *Bioinformatics* 23 (1): 127–28. doi:[10.1093/bioinformatics/btl529](https://doi.org/10.1093/bioinformatics/btl529).

Paradis, Emmanuel, Julien Claude, and Korbinian Strimmer. 2004. “APE: Analyses of Phylogenetics and Evolution in R Language.” *Bioinformatics* 20 (2): 289–90. doi:[10.1093/bioinformatics/btg412](https://doi.org/10.1093/bioinformatics/btg412).

Pennell, Matthew W, Jonathan M Eastman, Graham J Slater, Joseph W Brown, Josef C Uyeda, Richard G FitzJohn, Michael E Alfaro, and Luke J Harmon. 2014. “Geiger V2.0: An Expanded Suite of Methods for Fitting Macroevolutionary Models to Phylogenetic Trees.” *Bioinformatics* 30 (15): 2216–8. doi:[10.1093/bioinformatics/btu181](https://doi.org/10.1093/bioinformatics/btu181).

Revell, Liam J. 2012. “Phytools: An R Package for Phylogenetic Comparative Biology (and Other Things).” *Methods Ecol. Evol.* 3 (2). Blackwell Publishing Ltd: 217–23. doi:[10.1111/j.2041-210X.2011.00169.x](https://doi.org/10.1111/j.2041-210X.2011.00169.x).

Schliep, Klaus Peter. 2011. “Phangorn: Phylogenetic Analysis in R.” *Bioinformatics* 27 (4): 592–93. doi:[10.1093/bioinformatics/btq706](https://doi.org/10.1093/bioinformatics/btq706).
