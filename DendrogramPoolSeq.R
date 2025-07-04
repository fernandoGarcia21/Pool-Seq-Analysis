setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

#.libPaths(c('/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

#load libraries
library(dplyr)
library(adegenet) #package devoted to the multivariate analysis of genetic markers data. (df2genind)
library(NAM) #nested association mapping - Alencar Xavier (Gdist)
library(dendextend)
library(parallelly)

#libraries for bootstrapping
library(future)
plan(multisession)
library(shipunov)


#library(shipunov, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")
#library(furrr, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")
#library(NAM, lib.loc = "C:/Users/dgarciac/AppData/Local/Programs/R/R-4.3.0/library")
#library(dendextend, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")


###########################################################
#                Functions
###########################################################




###########################################################
#                Main body: Call to functions
###########################################################

#Input file paths

#Input file paths
frequencies_file = 'Datasets/subsetEurope.frq' #SNP-alleles frequencies file
sample_location_file = 'Datasets/Sample_Location.txt' #Two column file with the sample name and the geographical location
snp_sample_size = 100000

exclude_populations = c()
exclude_samples = c("17TFA",
                    "VA-3191-SPCAN-low_S22_L003",
                    "VA-3191-SPCAN-high_S23_L003",
                    "SPn_Crab_High",
                    "SPn_Wave_Low",
                    "SPs_Crab_High",
                    "SPs_Wave_Low",
                    "12POR")

#Read file and load data
frequencies_data_original = read.table(frequencies_file, header = TRUE, check.names = FALSE)
sample_location_data = read.table(sample_location_file, header = TRUE, sep = "\t", check.names = FALSE)

col_samples_start = 4 #The column number of the first sample
n_bootstrap_replicates = 200 #Number of dendrogram replicates to estimate bootstrap

#Find the geographic location according to the sample name, to exclude unwanted populations and unwanted samples
aux_n_cols = dim(frequencies_data_original)[2]
sample_names = data.frame(Sample = colnames(frequencies_data_original)[col_samples_start:aux_n_cols])
sample_names$Location = sample_location_data$LOCATION[match(sample_names$Sample, sample_location_data$SAMPLE)]

keep_samples = sample_names[!sample_names$Sample %in% exclude_samples, 'Sample'] #Find the rows (location names) that are not in the exclude list
genotypes_data = frequencies_data_original[,c(keep_samples)] #Filter the frequencies data to exclude unwanted locations 
rownames(genotypes_data) = paste(frequencies_data_original$chrom, frequencies_data_original$pos, sep = '_')

#Transpose the frequencies matrix to have populations in the rows and SNPs in columns
t_frequencies_data = t(genotypes_data)

#Subsample the columns to not use all SNPs
n_t_cols = ncol(t_frequencies_data)
tmp_random_snps = sample(n_t_cols, snp_sample_size, replace = FALSE)
s_genotypes_data = genotypes_data[, tmp_random_snps]

print("Computing genetic distances")
#Compute genetic distance using NAM package
gdist1 = Gdist(s_genotypes_data, method = 1) #This function computes measures of genetic distances between populations using a genpop object. 

#Create multiple dendrograms from random sampling with replacement, for bootstrapping.
#The process occurs in para llel thanks to the furrr::future_map function
list_of_hc <- local({
  options(future.globals.maxSize = 1024 * 1024 * 1024) # Adjust the size as needed (e.g., 1GB)
  
  furrr::future_map(1:n_bootstrap_replicates, function(i) {
    ##create a dataframe with replacement using original df, subsampling 100% of the rows
    tmp_subset = as.data.frame(s_genotypes_data) %>% slice_sample(prop = 1, replace=TRUE)
    
    ##run hclust on the data
    dist_mat = NAM::Gdist(as.matrix(tmp_subset), method = 1) # Assuming Gdist is from the NAM package.
    hc = hclust(dist_mat, method = 'ward.D') #Performs hierarchical clustering on the distance matrix of the sample
    ##save the hclust result to a list
    hc
  }, .progress = TRUE,
  .options = furrr::furrr_options(packages = "NAM", seed = TRUE))
})


#first element of the list is based on original df
list_of_hc[[1]] <- hclust(gdist1,method='ward.D')


print("Plotting the dendrogram")
# Create dendrogram
#Create a tree plot object to extract
bb3 <- Bclust(hclist=list_of_hc, relative = TRUE)
plt_bclust = plot(bb3) #Plot object from Bclust that we will use to extract the bootstrap coordinates

#Transform the tree into a dendrogram
dend = as.dendrogram(bb3$hclust)


# Change the size of the labels
dend <- set(dend, "labels_cex", 0.5)

locations_dend = sample_names$Location[match(labels(dend), sample_names$Sample)]

labels_colors(dend) <- as.numeric(as.factor(locations_dend))

location_colors = c("#8BBABB", "#a051b1", "#81B214", "#99820d", 
                    "#F7418F", "#9078cb", "#c5b61f", "#303997", 
                    "#449860", "#7754c0", "#95efe1", "#e9abff", 
                    "#FFA900", "#7c8d9a", "#a01c54", "#127f35", 
                    "#af281e", "#af5e00")

palette(location_colors)
print("Plotting figure")
#pdf('Dendrogram_Europe.pdf',width = 10, height = 6)
# plot dendrogram
plot(dend, horiz=T)
#Extract the bootstrapt coordinates from the Bclust plot and plot the labels on the dendrogram
text(x=plt_bclust$coords[,'y']+0.009, y=plt_bclust$coords[,'x']+0.85, 
     labels = ceiling(bb3$values*100), cex = 0.5, col = "black")  # Plot bootstrap values
legend("topleft", legend=unique(as.factor(locations_dend)), 
       pch=16, 
       horiz = FALSE,
       cex=0.6,
       ncol = 2,
       col=unique(as.factor(locations_dend)))

#dev.off()

print("Dendrogram generated successfully and exported to Dendrogram_Europe.pdf")


#Write a file with the genetic distance between the samples and the average Venice
write.table(as.matrix(gdist1), file = 'scripts/Distance_matrix.txt', quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)



