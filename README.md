**DendrogramPoolSeq.R:** This script is used to generate dendrograms from Pool-Seq data, likely for visualizing genetic relationships between samples or populations. It utilizes libraries such as dendextend, adegenet, and NAM for hierarchical clustering and plotting.

**DiversityDistributionsPoolSeq.R:** This script analyzes and visualizes diversity statistics (like FST, Pi, and Tajima's D) across different populations. It includes functions for reading and merging data files and generating various plots, potentially including distributions of these statistics.

**DriftSimulation.R:** This script performs simulations related to genetic drift and population size changes over generations. It includes a function (population_logistic_growth) to estimate population size based on logistic growth in a neutral demographic inference model.

**DxyAndFSTGerendalfPoolSeq.R:** This script focuses on calculating and visualizing Dxy (absolute genetic divergence) and FST (fixation index) between populations. It generates scatter plots of FST versus Dxy, highlighting potential outlier regions.

**ExploreInversionsPoolSeq.R:** This script is designed to explore and visualize genomic inversions by plotting allele frequencies between populations. It includes functions to generate unique pairs of items (populations) and create scatter plots to identify diagnostic SNPs within inversions.

**FindOutliersFSTAndAbasolutePi.R:** This script identifies genomic regions that are outliers for FST and absolute Pi (nucleotide diversity). It includes functions to generate Manhattan plots for these statistics and likely uses permutation tests or other statistical methods to define outliers.

**GenomeWideFSTPoolSeq.R:** This script computes and visualizes genome-wide FST values between two populations. It contains a compute_fst function and generates Manhattan plots to display FST values across the genome.

**GOEnrichmentTopGo.R:** This script performs Gene Ontology (GO) enrichment analysis using the topGO package. It reads gene annotation files (GAF) and gene lists to identify enriched GO terms in specific gene sets, likely from outlier regions or genes under selection.

**GWS_FST_GerendalfPoolSeq.R:** This script generates genome-wide Manhattan plots for FST and Dxy, similar to GenomeWideFSTPoolSeq.R and DxyAndFSTGerendalfPoolSeq.R but specifically structured to combine and plot results for different population comparisons.

**GWS_Pi_GerendalfPoolSeq.R:** This script generates genome-wide Manhattan plots for nucleotide diversity (Theta Pi). It reads and merges diversity files and visualizes the distribution of Pi across chromosomes for different populations.
