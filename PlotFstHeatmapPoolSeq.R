setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

#install.packages('pheatmap')
#install.packages("pheatmap", lib = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3/")
library("pheatmap")
library(ggplot2)
library(reshape2) # For melting the matrix
library(dplyr)
library(grid) # To explicitly control page breaks
library(RColorBrewer)

######################################################
# Computes the FST for each SNP of two populations 
# and returns an array with FST values
######################################################
compute_fst = function(aaf_pop_A, aaf_pop_B){
  raf_pop_A = 1 - aaf_pop_A
  raf_pop_B = 1 - aaf_pop_B
  
  expec_heterozygosity_A = 2 * raf_pop_A * aaf_pop_A
  expec_heterozygosity_B = 2 * raf_pop_B * aaf_pop_B
  
  h_s = (expec_heterozygosity_A + expec_heterozygosity_B) / 2
  h_t_ref = (raf_pop_A + raf_pop_B) / 2
  h_t_alt = (aaf_pop_A + aaf_pop_B) / 2
  h_t = 2 * h_t_ref * h_t_alt
  #Find the indexes of all ht equal 0
  h_t_0 = which(h_t == 0)
  
  fst = (h_t - h_s) / h_t
  #print(h_t_0)
  #Replace those Nan values of fst with 0
  fst[h_t_0] = 0
  
  return(fst)
}


#Estimate and write inversion frequencies for all populations
file_af_snps = 'Datasets/subsetEurope.frq';
output_file_name = 'Europe_FST_Heatmap.pdf';
col_start_freq = 4


args = commandArgs(trailingOnly=TRUE)
# test if there is at least five arguments: if not, return an error
if (length(args) < 2) {
  stop("At least two arguments must be supplied (allele frequencies file and output file name)", call.=FALSE)
} else {
  file_af_snps = args[1] #Path of the file where the AF of the pools is found
  output_file_name = args[2] #Path of the file where the af of the diagnostic SNPs
}

print('Finish reading parameters')


#Read the diagnostic SNPs af
data_af_snps = read.table(file_af_snps, header = TRUE, check.names = FALSE)
#Preserve the meta data of the SNPs
data_SNP_metainfo = data_af_snps[,1:col_start_freq-1]
#Remove the meta data from the AF dataframe to keep only the frequencies columns
data_af_snps = data_af_snps[,col_start_freq:dim(data_af_snps)[2]]

# 1. Get the names of the populations
population_names = colnames(data_af_snps)
num_populations = length(population_names)

# 2. Generate all unique pairs of populations
population_pairs = combn(population_names, 2)

# 3. Initialize a matrix or data frame to store the pairwise FST values
fst_matrix = matrix(NA, nrow = num_populations, ncol = num_populations)
colnames(fst_matrix) = population_names
rownames(fst_matrix) = population_names

# Alternatively, to store as a data frame:
fst_df = data.frame(Pop1 = character(), Pop2 = character(), FST = numeric(), stringsAsFactors = FALSE)

# 4. Iterate through the population pairs and calculate FST
for (i in 1:ncol(population_pairs)) {
  pop1_name = population_pairs[1, i]
  pop2_name = population_pairs[2, i]
  
  # Extract allele frequencies for the two populations
  pop1_freq = data_af_snps[[pop1_name]]
  pop2_freq = data_af_snps[[pop2_name]]
  
  # Calculate FST using your function
  fst_values = compute_fst(pop1_freq, pop2_freq)
  mean_fst = mean(fst_values)
  
  # Store the result in the matrix
  fst_matrix[pop1_name, pop2_name] = mean_fst
  fst_matrix[pop2_name, pop1_name] = mean_fst # Matrix is symmetric
  
  # Store the result in the data frame
  fst_df = rbind(fst_df, data.frame(Pop1 = pop1_name, Pop2 = pop2_name, FST = mean_fst))
}

#################################################################
#Plot the heatmap
#################################################################

# 1. Convert the FST matrix to a long format (tidy data)
fst_melted = melt(fst_matrix, na.rm = TRUE)
colnames(fst_melted) = c("Population1", "Population2", "FST")


# Determine the order of populations based on FST
# We need to find a way to order the populations based on their overall FST values
# in the lower triangle. One approach is to calculate the average FST
# for each population across the comparisons in the lower triangle.

# Get unique population names
all_populations <- unique(c(as.character(fst_melted$Population1), as.character(fst_melted$Population2)))

# Calculate the average FST for each population in the lower triangle
population_fst_summary <- fst_melted %>%
  group_by(Population1) %>%
  summarise(mean_fst = mean(FST)) %>%
  rename(Population = Population1) %>%
  bind_rows(
    fst_melted %>%
      group_by(Population2) %>%
      summarise(mean_fst = mean(FST)) %>%
      rename(Population = Population2)
  ) %>%
  group_by(Population) %>%
  summarise(overall_mean_fst = mean(mean_fst)) %>%
  arrange(overall_mean_fst, decreasing = TRUE)

# Get the ordered list of populations
ordered_populations <- population_fst_summary$Population

# Convert Population1 and Population2 to ordered factors
fst_melted$Population1 <- factor(fst_melted$Population1, levels = ordered_populations)
fst_melted$Population2 <- factor(fst_melted$Population2, levels = ordered_populations)


# Create the heatmap using ggplot2 with sorted axes
plt_heatmap = ggplot(fst_melted, aes(x = Population1, y = Population2, fill = FST)) +
  geom_tile(color = "white") +
  #scale_fill_viridis_c(option = "plasma", name = "FST") +
  scale_fill_gradientn(colors = brewer.pal(n = 11, name = "BrBG"), name = "FST") + # Use scale_fill_gradientn
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "right") +
  coord_equal()



#Create a heatmap that includes a dendrogram using the euclidean distances of the FST estimats
plt_heatmapTree = pheatmap(fst_matrix, cutree_rows = 4, cutree_cols = 4, color = brewer.pal(n = 11, name = "BrBG"))




pdf(output_file_name, width = 10, height = 8)
print(plt_heatmap)
# Force a new page before plotting pheatmap
grid.newpage()
print(plt_heatmapTree)



#Create a barplot with the FST of Venice or San Francisco vs the other populations
focus_populations = c('VEN','RED','OAK')
fst_focus_pops = fst_df[fst_df$Pop1 %in% focus_populations, c(1,2,3)]
tmp_second_part = fst_df[fst_df$Pop2 %in% focus_populations, c(2,1,3)]
names(tmp_second_part) = names(fst_focus_pops)
fst_focus_pops = rbind(fst_focus_pops, tmp_second_part)
unique_focus_pops = unique(fst_focus_pops$Pop1)


for(tmp_pop in unique_focus_pops){
  
  fst_focus_pop = fst_focus_pops[fst_focus_pops$Pop1 == tmp_pop,]
  fst_focus_pop = fst_focus_pop[order(fst_focus_pop$FST, decreasing = TRUE), ]
  
  fst_focus_pop$Pop2 <- factor(fst_focus_pop$Pop2, levels = unique(fst_focus_pop$Pop2))
  
  plt_barplot = ggplot(fst_focus_pop, aes(x = Pop2, y = FST, fill = FST)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = brewer.pal(n = 11, name = "BrBG"), name = "FST") + # Using a fixed number of colors for the gradient
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(), # Hide the x-axis label
          plot.title = element_text(hjust = 0.5) # Center the title (optional)
    )+
    labs(title = tmp_pop) # Add the plot title
  
  # Force a new page before plotting the barplot
  grid.newpage()
  print(plt_barplot)
}

dev.off()

print('Output file created')





