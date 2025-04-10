All codes uploaded here are associated with the manuscript "Selecting genes for phylogenetic analysis based on the geometry of gene-treespace" by Xiumei Lu et al.

property_calculation: A document containing scripts for calculating gene properties.

1_tree_dist_tutorial_Quartet_2041.r: A script for calculating the quartet distance matrix and performing NMDS (supplementary Fig. 2).

2_property_space.r: A script used for calculating PC2 based on the gene property matrix (supplementary Fig. 1).

3_bins.r: A script for comparing properties by mapping selected gene trees on treespace (Fig. 3; supplementary Figs. 3-17).

4_why1000.r: A script for determining the selection of 1,000 genes (Supplementary Figs. 18-33).

5_occupation_1000.r: A script testing if properties of phylogenetically informative genes coincide with the center of treespace, using 1,000 gene subsamples (Fig. 3).

6_occupation_500.r: A script testing if properties of phylogenetically informative genes coincide with the center of treespace, using 500 gene subsamples (supplementary Fig. 34).

7_permutation_test_1000.r: A script generating permutation tests and comparing the number of shared trees in 1,000-gene sets to random sets (Fig. 4).

8_permutation_test_500.r: A script generating permutation tests and comparing the number of shared trees in 500-gene sets to random sets (supplementary Fig. 36).

9_subsets_concatentations_faa.sh: A script for concatenating genes selected by various criteria; the script catfasta2phyml.pl (attached) is cited.

10_facetoface_tree_visual.r: A script for comparing phylogeny from the center with various criteria (Figs. 6-7; supplementary Figs. 36-53).

11_TreeSimiCalculandRank.r: A script calculating the similarity of trees from various criteria to those from the center, and plotting ranks (supplementary Fig. 54).

12_tree_dist_tutorial_Quartet_20.r: A script for calculating the quartet distance matrix of trees from various criteria and the center.

13_HierarchicalClusteringHeatmap.py: A script for generating a hierarchical cluster map (Fig. 8).

14_AU_test.sh: A script for conducting the AU test.

15_tree_dist_tutorial_Quartet_2062.r: A script for calculating the quartet distance matrix and performing PCA.

16_TreespacewithcounterlineZoomin.r: A script for plotting Fig. 1 using quartet distance from the above script.

Functions.R: A script containing functions for generating quartet distance matrix.

catfasta2phyml.pl: A script for gene concatenation.
