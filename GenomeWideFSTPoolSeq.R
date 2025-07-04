setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)


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
focus_population = 'VEN'
list_pops_vs = c('3GAL','2FAR','1BER')
window_size = 10000 #Window size in base positions along the genome
output_file_name = 'Europe_FST_Manhattan.jpg';
col_start_freq = 4


args = commandArgs(trailingOnly=TRUE)
# test if there is at least five arguments: if not, return an error
if (length(args) < 5) {
  stop("At least two arguments must be supplied (allele frequencies file, focus population like VEN, and output file name)", call.=FALSE)
} else {
  file_af_snps = args[1] #Path of the file where the AF of the pools is found
  focus_population = args[2] #E.g. VEN
  aux_pops_ver = args[3] #Multiple populations separated by comma
  window_size = as.numeric(args[4]) #Sliding window size in basepairs
  output_file_name = args[5] #Path of the file where the af of the diagnostic SNPs
  
  list_pops_vs = strsplit(aux_pops_ver, ',')[[1]]
}

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




#Prepare the metadata to plot the X axis

#Extract the LG number from the LG name
aux_df_row_names$LG = sapply(strsplit(aux_df_row_names$chrom, "_"), function(x) x[1])
aux_df_row_names$LGN = as.numeric(substr(aux_df_row_names$LG,3,length(aux_df_row_names$LG)))

#Merge the FST dataset with the metadata
data_fst_plot = as.data.frame(cbind(fst_mains_matrix, aux_df_row_names))
rownames(data_fst_plot) = rownames(fst_mains_matrix)

#Sort the matrix by the LG number and window number
data_fst_plot = data_fst_plot[order(data_fst_plot$LGN, data_fst_plot$window),]


#Compute the cumulative of the windows to plot each value after the other
fst_df_sorted_final = data_fst_plot %>% 
  # Compute chromosome size
  group_by(LGN) %>% 
  summarise(LG_len=max(window)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(LG_len)-LG_len) %>%
  select(-LG_len) %>%
  
  # Add this info to the initial dataset
  left_join(data_fst_plot, ., by=c("LGN"="LGN")) %>%
  
  # Add a cumulative position of each SNP
  arrange(LGN, window) %>%
  mutate( BPcum=window+tot)


# Determine the coordinates of the shaded rectangles and the x position of the axis labels
axis_fst_df <- fst_df_sorted_final %>% 
  group_by(LGN) %>% 
  summarize(start=min(BPcum),
            end=max(BPcum))
#Add the gray white intercalated
axis_fst_df$color <- rep(c(1,0),ceiling(nrow(axis_fst_df)/2))[1:nrow(axis_fst_df)]

#Estimate the center of the chromosome
ref_end_previouslg = append(0, axis_fst_df$end+0.0001)
ref_end_previouslg = ref_end_previouslg[1:length(ref_end_previouslg)-1]
axis_fst_df$start = ref_end_previouslg
axis_fst_df$center = axis_fst_df$start + ((axis_fst_df$end - axis_fst_df$start) / 2)


#Create a list of Manhattan plots
list_manhattan_plt = list()
x_pos_percentiles = (axis_fst_df[axis_fst_df$LGN %in% c(5,6),])[1,]$start + window_size

#For loop iterating over the populations in the list_pops_vs
for(pop_vs in list_pops_vs){
  
  fst_df_pop = fst_df_sorted_final[,c('BPcum', pop_vs)]
  colnames(fst_df_pop) = c('pos_acum','fst')
  
  # Calculate the 95th and 99th percentile
  percentile_95 = quantile(fst_df_pop$fst, probs = 0.95)
  percentile_99 = quantile(fst_df_pop$fst, probs = 0.99)
  
  plt = ggplot()+
    # Show all points
    geom_point(data= fst_df_pop, aes(x=pos_acum, y=fst), size = 0.5, color="#1d3557", alpha=0.7) +
    # custom X axis:
    scale_x_continuous( label = axis_fst_df$LGN, breaks= axis_fst_df$center , expand = c(0.01,0)) +
    #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Add the horizontal red line at the 95th percentile
    geom_hline(yintercept = percentile_95, color = "blue", linetype = "dashed", linewidth = 0.8) +
    # Add a text label above the line
    geom_text(
      data = data.frame(
        x = x_pos_percentiles,
        y = round(percentile_95, 2),
        label_text = paste0(round(percentile_95, 2), " (95% percentile)")
      ),
      aes(x = x, y = y, label = label_text),
      hjust = -0.1, # Adjust horizontal justification
      vjust = -0.5, # Adjust vertical justification
      color = "blue",
      size = 4
    ) +
    # Add the horizontal red line at the 95th percentile
    geom_hline(yintercept = percentile_99, color = "red", linetype = "dashed", linewidth = 0.8) +
    # Add a text label above the line
    geom_text(
      data = data.frame(
        x = x_pos_percentiles,
        y = round(percentile_99, 2),
        label_text = paste0(round(percentile_99, 2), " (99% percentile)")
      ),
      aes(x = x, y = y, label = label_text),
      hjust = -0.1, # Adjust horizontal justification
      vjust = -0.5, # Adjust vertical justification
      color = "red",
      size = 4
    ) +
    
    
    # Custom the theme:
    theme_minimal() +
    theme(
      legend.position="none",
      legend.title = element_text(size = 8, face = "bold"),
      panel.border = element_blank(),
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
    theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"))+ # Top, Right, Bottom, Left
    #ylim(0.0, 1.0)+
    labs(y = "FST", x="") +
    ggtitle(paste(focus_population, 'vs', pop_vs))+
    geom_rect(data=axis_fst_df, mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=as.factor(color)), alpha=0.15, show.legend = FALSE)+
    scale_fill_manual(values = c('1'="gray", '0'="white"), name=NULL)+
    guides(fill = FALSE)
  
  list_manhattan_plt[[pop_vs]] = plt
  
}

#Generate a PDF with the Manhattan plots on different rows

n_rows_plot = length(list_pops_vs)

width_page = 15
height_page = 4 * n_rows_plot #Inches
format = file_ext(output_file_name)

if(format == 'pdf' | format == 'PDF'){
  pdf(output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_manhattan_plt, ncol = 1, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_manhattan_plt, ncol = 1, nrow = n_rows_plot, align = "v")
    
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
  


#Estimate nucleotide diversity pi as pi = 2pq https://onlinelibrary.wiley.com/doi/10.1111/mec.12522


  