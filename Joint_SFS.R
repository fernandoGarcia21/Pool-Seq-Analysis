setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')


#########################################
# Creates the joint site frequency spectrum for two given populations
# @author: Diego Garcia
# @date: 12/05/2025
#########################################

library(rlang) # For the !! operator
library(foreach)
library(doParallel)
library(parallelly)
library(ggplot2)

file_derived_allele_counts = 'Datasets/JSFS/counts.derivedAllele.collinear.SanFrancisco.sync'
file_sites_coverage = 'Datasets/JSFS/coverage.derivedAllele.collinear.SanFrancisco.sync'
file_populations_coverage = 'Datasets/JSFS/coverage.populations.txt'

tmp_pop1 = 'OAK'
tmp_pop2 = 'SM-BOO-high-merged'

#args = commandArgs(trailingOnly=TRUE)
# test if there is at least five arguments: if not, return an error
#if (length(args) < 5) {
#  stop("At least two arguments must be supplied", call.=FALSE)
#} else {
#  file_derived_allele_counts = args[1] #Path of the file where the allele counts of the pools is found
#  file_sites_coverage = args[2] #Path of the file with the coverage of the SNPs
#  file_populations_coverage = args[3] #Path of the file with the coverage of the populations
#  tmp_pop1 = args[4] #Name of the population that will appear on X axis
#  tmp_pop2 = args[5] #Name of the population that will appear on Y axis
#}

print('Finish reading parameters')
print(paste(tmp_pop1, 'vs', tmp_pop2))


#Read the input datasets
data_derived_allele_counts = read.table(file_derived_allele_counts, header = TRUE, sep = '\t', check.names = FALSE)
data_sites_coverage = read.table(file_sites_coverage, header = TRUE, sep = '\t', check.names = FALSE)
data_populations_coverage = read.table(file_populations_coverage, header = TRUE, sep = '\t', check.names = FALSE)

#Exclude SNPs that have larger allele counts than the average coverage of the population
#idx_exclude = which(data_derived_allele_counts[,tmp_pop1] > data_populations_coverage[data_populations_coverage$`SAMPLE NAME` == tmp_pop1, ]$Coverage |
#                      data_derived_allele_counts[,tmp_pop2] > data_populations_coverage[data_populations_coverage$`SAMPLE NAME` == tmp_pop2, ]$Coverage)
#data_derived_allele_counts = data_derived_allele_counts[-idx_exclude, ]


############################################
# For development, trim the dataset to use only a subset
data_derived_allele_counts = data_derived_allele_counts[data_derived_allele_counts$chrom %in% c('LG5_SUPER_9', 'LG12_SUPER_3'),]
data_sites_coverage = data_sites_coverage[data_sites_coverage$chrom %in% c('LG5_SUPER_9', 'LG12_SUPER_3'),]
############################################



#Normalize the the derived alleles based on the coverage by sampling

#Check whether the populations must be downsampled
tmp_coverage_pop1 = data_populations_coverage[data_populations_coverage$`SAMPLE NAME` == tmp_pop1, ]$Coverage
tmp_coverage_pop2 = data_populations_coverage[data_populations_coverage$`SAMPLE NAME` == tmp_pop2, ]$Coverage

downsampled_population = ''
target_coverage_population = ''

downsampling_required = FALSE
#When the population 1 has larger coverage, the SNPs must be downsampled to the population 2 coverage per snp
if(tmp_coverage_pop1 > tmp_coverage_pop2){
  downsampled_population = tmp_pop1
  target_coverage_population = tmp_pop2
  downsampling_required = TRUE
}else{
  if(tmp_coverage_pop2 > tmp_coverage_pop1){
    downsampled_population = tmp_pop2
    target_coverage_population = tmp_pop1
    downsampling_required = TRUE
  }
}

#Check if there is downsampling required otherwise just use the given allele counts
if(downsampling_required) {
  
  avg_pop_size = round(mean(c(tmp_coverage_pop1, tmp_coverage_pop2)))
  print(paste('The average population size will be', avg_pop_size))
  
  #Get the number of chromosomes and parallelize the normalization by chromosome
  unique_chromosomes = unique(data_derived_allele_counts$chrom)
  
  # Downsample the counts of the population with the largerst coverage
  # in parallel to make it more efficient
  cores=availableCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  #For each chromosome, downstream the counts and merge the results as a dataframe
  downsampled_counts <- foreach(
    chr = unique_chromosomes,
    .combine = 'rbind'
  ) %dopar% {
    
    #Object that will store the dataset after downsampling
    downsampling_counts = data_derived_allele_counts[data_derived_allele_counts$chrom == chr, c('chrom','pos',tmp_pop1, tmp_pop2)]
    downsampling_coverage = data_sites_coverage[data_sites_coverage$chrom == chr, c('chrom','pos',tmp_pop1, tmp_pop2)]
    
    for (i_pos in downsampling_counts$pos) {
      
      # Estimate the frequency of the derived allele, it will be used as the downsampling probability
      downsample_prob_pop1 = downsampling_counts[downsampling_counts$pos == i_pos, tmp_pop1] / downsampling_coverage[downsampling_coverage$pos == i_pos, tmp_pop1]
      downsample_prob_pop2 = downsampling_counts[downsampling_counts$pos == i_pos, tmp_pop2] / downsampling_coverage[downsampling_coverage$pos == i_pos, tmp_pop2]
      
      # Downsample allele counts for both populations using binomial distribution
      downsampling_counts[downsampling_counts$pos == i_pos,tmp_pop1] = rbinom(1, avg_pop_size, downsample_prob_pop1)
      downsampling_counts[downsampling_counts$pos == i_pos,tmp_pop2] = rbinom(1, avg_pop_size, downsample_prob_pop2)
    }
    
    downsampling_counts
    
  }
  stopCluster(cl) #Stop the cores cluster
  
  
}else{
  print('Downsampling is not required')
  downsampled_counts = data_derived_allele_counts[, c('chrom','pos',tmp_pop1, tmp_pop2)]
}


#Create a dataframe with the summary of the pairwise allele counts frequencies
jsfs_plot_df = as.data.frame(table(downsampled_counts[,c(tmp_pop1, tmp_pop2)]))
colnames(jsfs_plot_df) = c(tmp_pop1, tmp_pop2, "Count")
jsfs_plot_df[[tmp_pop1]] <- as.numeric(as.character(jsfs_plot_df[[tmp_pop1]]))
jsfs_plot_df[[tmp_pop2]] <- as.numeric(as.character(jsfs_plot_df[[tmp_pop2]]))

jsfs_plot_df_filtered = jsfs_plot_df[jsfs_plot_df$Count >0,]
jsfs_plot_df_filtered = jsfs_plot_df_filtered[2:nrow(jsfs_plot_df_filtered),] #Remove the bottom left corner to see more colours
#jsfs_plot_df_filtered = jsfs_plot_df_filtered[jsfs_plot_df_filtered$OAK < 80 & jsfs_plot_df_filtered$`SM-BOO-high-merged` < 80, ]
jsfs_plot_df_filtered$LogCount = log(jsfs_plot_df_filtered$Count, base = 10) #Scale the allele number for visualization


# Plot the JSFS using ggplot2
 ggplot(jsfs_plot_df_filtered, aes(x = !!sym(tmp_pop1), y = !!sym(tmp_pop2), fill = LogCount)) +
            geom_tile() +
            scale_fill_gradientn(colours = rainbow(20)) + # Adjust the number in rainbow() for more or fewer colors
            labs(x = paste("Derived Allele Count in ", tmp_pop1),
                 y = paste("Derived Allele Count in ", tmp_pop2),
                 fill = "Number of SNPs",
                 title = "Joint Site Frequency Spectrum")+
            #xlim(0,200)+
            #ylim(0,200)+
            theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),)

###########################################################
#Print the plot as an image
output_file_name = paste('JSFS',tmp_pop1, 'vs', tmp_pop2,'jpg', sep = '.')
ggsave(output_file_name, 
       plot = plt_jsfs, 
       device = "jpeg", 
       units = "in",
       width = 9, 
       height = 8, 
       dpi = 300)
###########################################################