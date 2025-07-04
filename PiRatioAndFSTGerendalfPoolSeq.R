setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
library(RColorBrewer)



read_and_merge_files = function(p_file_list, pool_names_array, p_sep){
  
  # Initialize an empty list to store the data frames
  data_frames_list <- list()
  col_names_statistics = c()
  # Loop through each file in the list
  for (file_path in p_file_list) {
    # Read the CSV file into a data frame
    current_df = read.table(file_path, header = TRUE, sep = p_sep, check.names = FALSE)
    column_names_df = colnames(current_df)
    index_cols_vs = unlist(lapply(pool_names_array, function(pattern) which(grepl(pattern, column_names_df, fixed = TRUE))))
    
    col_names_statistics = column_names_df[index_cols_vs]
    current_df = current_df[,c(1,2,3, index_cols_vs)]
    
    # Add the current data frame to the list
    data_frames_list[[length(data_frames_list) + 1]] = current_df
    
    
    # Optional: Print the name of the file that was just loaded
    cat("Loaded file:", basename(file_path), "\n")
  }
  
  # Merge all data frames in the list into a single data frame
  if (length(data_frames_list) > 0) {
    merged_df = do.call(rbind, data_frames_list)
    cat("\nSuccessfully merged", length(data_frames_list), "files into a single data frame.\n")
    
    colnames(merged_df) = c("chrom", "start", "end", col_names_statistics)
  } else {
    cat("\nNo .fst-list.csv files found in the specified folder.\n")
    merged_df = NULL # Or however you want to handle the case with no files
  }
  return(merged_df)
}



#Directory where the diversity statistics are found
#fst_files_dir_path = 'Grenedalf/WindowFST Europe';
files_dir_fst_path = 'Grenedalf/WindowFST_20K';
files_dir_diversity_path = 'Grenedalf/WindowDiversity_20K';

#output_file_name = 'Europe_PiRatioFST_Grenedalf.pdf';
output_file_name = 'NonBiased.PiRatioFST_Grenedalf.jpg';

list_pops_vs = c('SM-BOO-low-merged',
                 'SM-BOO-high-merged',
                 'PBY_BLMP1',
                 'PTY_BLMP2',
                 'PLY-high'
                 )
#list_pops_vs = c('EAS-low',
#                 'T_NC_NO-high',
#                 'Fr_Crab_Low',
#                 'EAS-high',
#                 '3GAL',
#                 '18VAD'
#)

focus_population = 'RED'
area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['VEN']] = 'EU'

y_offset = list()
y_offset[['RED']] = -11
y_offset[['OAK']] = -11
y_offset[['VEN']] = -7

output_file_name = paste(focus_population, output_file_name, sep = '.')

# List all files in the specified folder with the .diversity.csv extension
file_list_fst = list.files(path = files_dir_fst_path, pattern = paste0("\\.",focus_population,".*\\.fst\\.csv$"), full.names = TRUE)
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)

merged_fst_df = read_and_merge_files(file_list_fst, list_pops_vs, '\t')
merged_diversity_df = read_and_merge_files(file_list_diversity, c(focus_population, list_pops_vs), '\t')

# Merge the dataframes by the first three columns
merged_df = merge(merged_fst_df, merged_diversity_df, by = c("chrom", "start", "end"))

print('Finish reading datasets')
print(paste('Snps in the FST dataset:', dim(merged_fst_df)[1]))
print(paste('Snps in the Diversity dataset:', dim(merged_diversity_df)[1]))
print(paste('Snps in the merged dataset (these SNPs will be used for the plot):', dim(merged_df)[1]))

#Delete the initial dataframes to release memmory
rm(merged_fst_df, merged_diversity_df)


#For loop iterating over the populations in the list_pops_vs
list_fst_pi_plt = list()
tmp_colnames_df = colnames(merged_df)
for(pop_vs in list_pops_vs){
  
  #Subset the dataframe to get only the fst and pi statistics for the current vs population
  tmp_col_fst = paste0( focus_population, ':', pop_vs, '.fst')
  tmp_col_pi_focus = paste( focus_population, 'theta_pi' , sep = '.')
  tmp_col_pi_vs = paste( pop_vs, 'theta_pi' , sep = '.')
  tmp_stats_pop_vs = merged_df[,c('chrom','start','end',tmp_col_fst,tmp_col_pi_focus,tmp_col_pi_vs)]
  tmp_stats_pop_vs = na.omit(tmp_stats_pop_vs) #Remove rows with NA values 
  colnames(tmp_stats_pop_vs) =  c('chrom','start','end','fst','pi_focus','pi_vs')
  
  #Estimate the pi ratio between the vs and the focus population, 
  # so that ratio > 1 indicates lower diversity in the focus population
  tmp_stats_pop_vs$pi_ratio = log(tmp_stats_pop_vs$pi_vs / tmp_stats_pop_vs$pi_focus)
  
  #Remove nan and infinite values (e.g. caused by pi = 0 in the focus population)
  tmp_stats_pop_vs = tmp_stats_pop_vs[!(is.nan(tmp_stats_pop_vs$pi_ratio) | is.infinite(tmp_stats_pop_vs$pi_ratio)), ]
 
  # Calculate the 95th percentile of the fst
  tmp_fst_percentile = 99
  tmp_pi_percentile = tmp_fst_percentile
  if(pop_vs == 'PLY-high'){
    tmp_pi_percentile = 99
  }
  percentile_fst = quantile(tmp_stats_pop_vs$fst, probs = tmp_fst_percentile/100, na.rm = TRUE)
  # Calculate the 95th and 5th percentile of the pi ratio
  percentile_pi_upper = quantile(tmp_stats_pop_vs$pi_ratio, probs = tmp_pi_percentile/100, na.rm = TRUE)
  percentile_pi_lower = quantile(tmp_stats_pop_vs$pi_ratio, probs = 1-tmp_pi_percentile/100, na.rm = TRUE)
  
  #Use the percentiles to assign a group (to colour in the plot) to the windows
  tmp_stats_pop_vs$group = 'Neutral'
  
  
  if(dim(tmp_stats_pop_vs[tmp_stats_pop_vs$fst >= percentile_fst & 
                          tmp_stats_pop_vs$pi_ratio >= percentile_pi_upper
                          , ])[1] > 0){
    
    tmp_stats_pop_vs[tmp_stats_pop_vs$fst >= percentile_fst & 
                       tmp_stats_pop_vs$pi_ratio >= percentile_pi_upper
                     , ]$group = 'Right'
  }
  
#  if(dim(tmp_stats_pop_vs[tmp_stats_pop_vs$fst >= percentile_fst & 
  #                          tmp_stats_pop_vs$pi_ratio <= percentile_pi_lower
  #                          , ])[1] > 0){
                            #    tmp_stats_pop_vs[tmp_stats_pop_vs$fst >= percentile_fst & 
                            #                       tmp_stats_pop_vs$pi_ratio <= percentile_pi_lower
                            #                     , ]$group = 'Left'
                            # }
  
  #print(dim(tmp_stats_pop_vs[tmp_stats_pop_vs$pi_ratio >= percentile_pi_upper, ])[1])
  #print(dim(tmp_stats_pop_vs[tmp_stats_pop_vs$pi_ratio <= percentile_pi_lower, ])[1])
  print('------------------------------------')
  
  max_pi = ceiling(mean(tmp_stats_pop_vs$pi_ratio) + sd(tmp_stats_pop_vs$pi_ratio)*10)
  min_pi = floor(min(tmp_stats_pop_vs$pi_ratio))
  print(max_pi)
  
  #Export the outlier windows
  tmp_outlier_file_name = paste(focus_population, pop_vs, 'WindowOutliers.NonBiased', 'txt', sep = '.')
  write.table(x = tmp_stats_pop_vs[tmp_stats_pop_vs$group == 'Right', ], file = tmp_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  #Generate a plot
  plt = ggplot()+
    # Show all points
    geom_point(data=tmp_stats_pop_vs, aes(x=pi_ratio, y=fst, col=group, alpha=0.7), size = 1.5, shape=19) +
    # Add the horizontal red line at the 95th percentile of the fst
    geom_hline(yintercept = percentile_fst, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
    # Add the vertical lines at the 95th and 5th percentiles
    geom_vline(xintercept = percentile_pi_upper, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
    #geom_vline(xintercept = percentile_pi_lower, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
    
    geom_text(
      data = data.frame(
        x = percentile_pi_upper,
        y = percentile_fst,
        label_text = paste0(nrow(tmp_stats_pop_vs[tmp_stats_pop_vs$group == 'Right', ]), " outlier windows")
      ),
      aes(x = x, y = y, label = label_text),
      hjust = 1.5, # Adjust horizontal justification
      vjust = y_offset[[focus_population]], # Adjust vertical justification
      color = "#80cdc1",
      size = 4.5
    ) +
    
    #ylim(-0.1, 0.6)+
    xlim(min_pi, max_pi)+
    # Custom the theme:
    theme_minimal() +
    theme(
      legend.position="none",
      legend.title = element_text(size = 8, face = "bold"),
      panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
      panel.grid.major.x = element_line(linewidth = 0.1),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.1),
      plot.title = element_text(hjust = 0.5, face='bold', size=12,),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=8),
      axis.text.x=element_text(size=8),
      plot.margin = margin( 1,0.0,0.0,0.0, "cm")
      
    )+
    scale_color_manual(values = c(Neutral = "#d4d4d4", Right = "#80cdc1", Left = "#a6611a"))+
    #scale_x_continuous(breaks = min_pi-1:max_pi)+ # Often the default shows all digits
    theme(plot.margin = unit(c(0.3, 0.5, 0.2, 0.5), "cm"))+ # Top, Right, Bottom, Left
    labs(y = "FST", x=paste0("log(PI ratio) (", pop_vs, "/", focus_population, ")")) +
    ggtitle(paste(focus_population, 'vs', pop_vs))
  
  
  list_fst_pi_plt[[pop_vs]] = plt
  
}



#Generate a PDF with the Scatter plots plots

n_cols_plot = 1
if(length(list_pops_vs) > 1){
  n_cols_plot = 3
}

n_rows_plot = ceiling(length(list_pops_vs)/n_cols_plot)

width_page = n_cols_plot * 5
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


