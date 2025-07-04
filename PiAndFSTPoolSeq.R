setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
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

# Function to apply (e.g., calculate the mean)
compute_pi <- function(x) {
  return(2*x*(1-x))
}


#Estimate and write inversion frequencies for all populations
file_af_snps = 'Datasets/subsetEurope.frq';
focus_population = 'VEN'
list_pops_vs = c('3GAL','2FAR','1BER')
window_size = 10000 #Window size in base positions along the genome
output_file_name = 'Europe_FST_PI_Scatter.jpg';
col_start_freq = 4


#args = commandArgs(trailingOnly=TRUE)
## test if there is at least five arguments: if not, return an error
#if (length(args) < 5) {
#  stop("At least two arguments must be supplied (allele frequencies file, focus population like VEN, and output file name)", call.=FALSE)
#} else {
#  file_af_snps = args[1] #Path of the file where the AF of the pools is found
#  focus_population = args[2] #E.g. VEN
#  aux_pops_ver = args[3] #Multiple populations separated by comma
#  window_size = as.numeric(args[4]) #Sliding window size in basepairs
#  output_file_name = args[5] #Path of the file where the af of the diagnostic SNPs
#  
#  list_pops_vs = strsplit(aux_pops_ver, ',')[[1]]
#}

print('Finish reading parameters')


#Read the diagnostic SNPs af
data_af_snps = read.table(file_af_snps, header = TRUE, check.names = FALSE)
#Preserve the meta data of the SNPs
data_SNP_metainfo = data_af_snps[,1:col_start_freq-1]
#Remove the meta data from the AF dataframe to keep only the frequencies columns
data_af_snps = data_af_snps[,c(focus_population, list_pops_vs)]


#Determine the window of the SNP according to the sliding windows size
data_SNP_metainfo$window = ceiling(data_SNP_metainfo$pos / window_size)
unique_windows = unique(data_SNP_metainfo[,c(1,4)])


# 1. Get the names of the populations
population_names = colnames(data_af_snps)
num_populations = length(population_names)

# 2. Generate all unique pairs of populations
all_other_populations = population_names[population_names != focus_population]
population_pairs = data.frame(Pop1 = rep(focus_population,length(all_other_populations)), Pop2 = all_other_populations)


# 3. Initialize a matrix or data frame to store the pairwise FST values
fst_mains_matrix = matrix(NA, nrow = nrow(unique_windows), ncol = num_populations-1)
colnames(fst_mains_matrix) = all_other_populations

aux_df_row_names = data_SNP_metainfo %>% group_by(chrom, window) %>%
  summarise(window_name = max(chrom))

rownames(fst_mains_matrix) = paste(aux_df_row_names$window_name , aux_df_row_names$window,sep = '_')

# 4. Iterate through the population pairs and calculate FST for each SNP
for (i in 1:nrow(population_pairs)) {
  
  pop1_name = population_pairs[i, 1]
  pop2_name = population_pairs[i, 2]
  
  # Extract allele frequencies for the two populations
  pop1_freq = data_af_snps[[pop1_name]]
  pop2_freq = data_af_snps[[pop2_name]]
  
  # Calculate FST using your function
  fst_values = compute_fst(pop1_freq, pop2_freq)
  
  #Calculate the mean FST within each window
  tmp_meta_fst = cbind(data_SNP_metainfo, fst_values)
  tmp_mean_fst_df = tmp_meta_fst %>% group_by(chrom, window) %>%
    summarise(mean_fst = mean(fst_values))
  
  # Store the result in the matrix
  fst_mains_matrix[, pop2_name] = tmp_mean_fst_df$mean_fst
  
}


#5. Estimate nucleotide diversity pi as pi = 2pq https://onlinelibrary.wiley.com/doi/10.1111/mec.12522

# Calculate the pi nucleotide diversity only for the focus population (e.g. VEN)
pi_means_array = c()

pi_values = compute_pi(data_af_snps[[focus_population]])
tmp_meta_pi = cbind(data_SNP_metainfo, pi_values)
tmp_mean_pi_df = tmp_meta_pi %>% group_by(chrom, window) %>%
  summarise(mean_pi = mean(pi_values))

#Append the mean nucleotide diversity to the matrix of the mean fst
fst_mains_matrix = cbind(fst_mains_matrix, mean_pi = tmp_mean_pi_df$mean_pi)


#Merge the FST and Pi dataset with the metadata
data_fst_pi_plot = as.data.frame(cbind(fst_mains_matrix, aux_df_row_names))


# 6. Identify FST Windows beyond the 99th percentile and plot the nucleotide diversity vs FST

#Create a list of scatter plots
list_fst_pi_plt = list()

#For loop iterating over the populations in the list_pops_vs
for(pop_vs in list_pops_vs){
  
  fst_df_pi_pop = data_fst_pi_plot[,c(pop_vs, 'mean_pi', 'chrom', 'window')]
  colnames(fst_df_pi_pop) = c('mean_fst', 'mean_pi', 'chrom', 'window')
  
  # Calculate the 99th percentile
  percentile_99 = quantile(fst_df_pi_pop$mean_fst, probs = 0.5)
  #subset the fst dataframe to get only the top percentile fst values
  tmp_perc_fst_pi_df = fst_df_pi_pop[fst_df_pi_pop$mean_fst >= percentile_99, ]
  

  plt = ggplot(data=tmp_perc_fst_pi_df, aes(x=mean_pi, y=mean_fst, fill=mean_fst, alpha=0.7))+
    # Show all points
    geom_point(size = 2.5, shape=21) +
    scale_fill_gradientn(colors = brewer.pal(n = 11, name = "BrBG"), name = "FST") + # Define the color gradient
    ylim(0, 1)+
    xlim(0, 0.5)+
    # Custom the theme:
    theme_minimal() +
    theme(
      legend.position="none",
      legend.title = element_text(size = 8, face = "bold"),
      panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5, face='bold', size=12,),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=8),
      axis.text.x=element_text(size=12),
      plot.margin = margin( 1,0.0,0.0,0.0, "cm")
      
    )+
    theme(plot.margin = unit(c(0.3, 0.5, 0.2, 0.5), "cm"))+ # Top, Right, Bottom, Left
    #ylim(0.0, 1.0)+
    labs(y = "FST", x="PI") +
    ggtitle(paste(focus_population, 'vs', pop_vs))
    
  
  list_fst_pi_plt[[pop_vs]] = plt
  
}



#Generate a PDF with the Scatter plots plots

n_cols_plot = 2
n_rows_plot = ceiling(length(list_pops_vs)/n_cols_plot)

width_page = 10
height_page = 4.5 * n_rows_plot #Inches
format = file_ext(output_file_name)

if(format == 'pdf' | format == 'PDF'){
  pdf(output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_fst_pi_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_fst_pi_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
    ggsave(output_file_name, 
           plot = combined_page, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')


