setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
library(RColorBrewer)
library(ggvenn)
library(VennDiagram)
library(tidyverse)

read_and_merge_files = function(p_file_list, p_sep){
  
  # Initialize an empty list to store the data frames
  data_frames_list = list()
  # Loop through each file in the list
  for (file_path in p_file_list) {
    # Read the CSV file into a data frame
    current_df = read.table(file_path, header = TRUE, sep = p_sep, check.names = FALSE)
    
    tmp_pop_vs = sapply(strsplit(file_path, "\\."), function(x) x[2])
    current_df$group = tmp_pop_vs
    
    
    # Add the current data frame to the list
    data_frames_list[[length(data_frames_list) + 1]] = current_df
    
    
    # Optional: Print the name of the file that was just loaded
    cat("Loaded file:", basename(file_path), "\n")
  }
  
  # Merge all data frames in the list into a single data frame
  if (length(data_frames_list) > 0) {
    merged_df = do.call(rbind, data_frames_list)
    cat("\nSuccessfully merged", length(data_frames_list), "files into a single data frame.\n")
  } else {
    cat("\nNo files found in the specified folder.\n")
    merged_df = NULL # Or however you want to handle the case with no files
  }
  return(merged_df)
}


#Directory where the ourliers of the FST and viersity diversity statistics are found
files_dir_shared_outliers_path = 'Grenedalf/WindowOutliers_FST_PI_10K';

#Directory where the FST windows are found, we just need this to know how many windows were used
files_dir_diversity_path = 'Grenedalf/WindowDiversity_10K';

output_file_name = 'TopLowPiDirection_Grenedalf10k.NonBiased.jpg';

vs_population = 'East'
focus_population = 'OAK'
list_pops_vs_area = list()
list_pops_vs_area[['Maine']] = c('SM-BOO-high-merged', 'SM-BOO-low-merged', 'PBY_BLMP1', 'PTY_BLMP2')
list_pops_vs_area[['East']] = c('SM-BOO-high-merged', 'SM-BOO-low-merged', 'PBY_BLMP1', 'PTY_BLMP2','PLY-high')
list_pops_vs_area[['PLY-high']] = c('PLY-high')
list_pops_vs_area[['SM-BOO-high-merged']] = c('SM-BOO-high-merged')
list_pops_vs_area[['SM-BOO-low-merged']] = c('SM-BOO-low-merged')
list_pops_vs_area[['PTY_BLMP2']] = c('PTY_BLMP2')

area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['PLY-high']] = 'NA'
area_focus_population[['Maine']] = 'NA'
area_focus_population[['East']] = 'NA'
area_focus_population[['SM-BOO-low-merged']] = 'NA'
area_focus_population[['VEN']] = 'EU'
area_focus_population[['EAS-low']] = 'EU'


output_file_name = paste(vs_population, focus_population, output_file_name, sep = '.')

#Read outliers file 
file_list_shared_outliers = list.files(path = files_dir_shared_outliers_path, pattern = paste0("^",vs_population,".\\", focus_population, "\\.LowPiThreshold.sharedFstLowPiRegionsOutlierWindows\\.txt$"), full.names = TRUE)
#file_list_shared_outliers = list.files(path = files_dir_shared_outliers_path, pattern = paste0("^",vs_population,".\\", focus_population, "\\.uniqueOutlierWindows\\.txt$"), full.names = TRUE)
if(length(file_list_shared_outliers) > 0){
  shared_outliers_df = read_and_merge_files(file_list_shared_outliers, '\t')
}else{
  shared_outliers_df = data.frame()
}

# List all files in the specified folder with the .diversity.csv extension
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)
merged_diversity_df = read_and_merge_files(file_list_diversity, '\t')
tmp_cols_pi = paste0( c(focus_population, list_pops_vs_area[[vs_population]]), '.theta_pi')
merged_diversity_df = merged_diversity_df[,c('chrom','start','end',tmp_cols_pi)]

# Merge the shared outliers and diversity dataframes by the first three columns
merged_outlier_windows_diversity = merge(shared_outliers_df, merged_diversity_df, by = c("chrom", "start", "end"))
colnames(merged_outlier_windows_diversity) = sapply(strsplit(colnames(merged_outlier_windows_diversity), "\\."), function(x) x[1])

#Create unique id for each window
merged_outlier_windows_diversity = merged_outlier_windows_diversity %>%
  mutate(unique_id = paste(chrom, start, end, sep = "_"))

#Re arrange the dataframe to eliminate unnecesary columns
merged_outlier_windows_diversity = merged_outlier_windows_diversity[, c(ncol(merged_outlier_windows_diversity), 6:ncol(merged_outlier_windows_diversity)-1)]


aux_max_pi = c()
list_plots = list()
is_first_plot = TRUE

colors_points = c('RED' = '#f30057',
                  'OAK' = '#f30057', 
                  'SM-BOO-high-merged' = '#80cdc1', 
                  'SM-BOO-low-merged' = '#80cdc1', 
                  'PLY-high' = '#80cdc1', 
                  'PBY_BLMP1' = '#80cdc1', 
                  'PTY_BLMP2' = '#80cdc1')



#Iterate over all windows to plot the diversity of the focus population vs the rest of the populations
for(i in 1:nrow(merged_outlier_windows_diversity)){
  
  tmp_window_df = merged_outlier_windows_diversity[i, ]
  
  # 1. Get the unique_id
  id = tmp_window_df$unique_id
  
  # 2. Extract RED value
  red_value = tmp_window_df[,focus_population]
  
  # 3. Select other columns for transformation (excluding unique_id and RED)
  other_cols = colnames(tmp_window_df)[!(names(tmp_window_df) %in% c("unique_id", focus_population))]
  
  # 4. Initialize an empty list to store results
  result_list = list()
  
  # 5. Iterate through each of the 'other_cols'
  for (col_name in other_cols) {
    # Create a temporary dataframe for the current pair
    temp_df_red = data.frame(
      unique_id = id,
      Pair = paste0(focus_population, '-', col_name),
      Value = red_value,
      Population = focus_population,
      stringsAsFactors = FALSE
    )
    
    temp_df_other = data.frame(
      unique_id = id,
      Pair = paste0(focus_population, '-', col_name),
      Value = tmp_window_df[,col_name],
      Population = col_name,
      stringsAsFactors = FALSE
    )
    
    # Combine RED and the current column's data
    result_list[[col_name]] = rbind(temp_df_other, temp_df_red)
  }
  
  # 6. Combine all temporary dataframes into a single result dataframe
  final_df = do.call(rbind, result_list)
  
  # 7. Reset row names (optional, but good practice)
  row.names(final_df) = NULL
  
  #Add a sequence of number to define the x axis
  #The labels will be the Population column
  final_df$x_pos = 1:nrow(final_df)
  
  # Create custom labels and their corresponding positions
  custom_x_breaks = final_df$x_pos
  custom_x_labels = final_df$Population
  
 
  aux_max_pi = append(aux_max_pi, max(final_df$Value))
  
  # Plotting with ggplot2
  tmp_plt = ggplot(final_df, aes(x = as.factor(x_pos), y = Value, color = Population, group = Pair)) +
    geom_line(linewidth = 0.5, color = 'black') +          # Add lines connecting the points within each 'Pair'
    geom_point(size = 1.5) + # Add points
    labs(
      title = merged_outlier_windows_diversity[i, ]$unique_id,
      x = element_blank(),
      y = "Theta pi",
      color = "Population"
    ) +
    scale_color_manual(values = colors_points)+
    scale_y_continuous(n.breaks = 6)+
    scale_x_discrete(
      breaks = custom_x_breaks, # Use the names of your custom_breaks_labels vector
      labels = custom_x_labels,        # Use the values of your custom_breaks_labels vector
      drop = FALSE                          # Keep all levels, even if not explicitly in breaks (optional)
    )+
    #ylim(0, 0.081)+ #Oakland
    ylim(0, 0.085)+ #Redwood
    theme_minimal() + # A clean theme
    theme(axis.text.x = if(is_first_plot) element_text(angle = 45, hjust = 1, size = 4) else element_text(angle = 45, hjust = 1, size = 4, color = 'white'),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 4),
          axis.title.y = if(is_first_plot) element_text(size = 4) else element_blank(),
          legend.position = 'none', 
          panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 5)
          ) # Rotate x-axis labels for readability
  
  list_plots[[merged_outlier_windows_diversity[i, ]$unique_id]] = tmp_plt
  is_first_plot = FALSE
}
print(max(aux_max_pi))
#Generate a PDF with the Scatter plots plots

n_cols_plot = 6
n_rows_plot = ceiling(length(list_plots)/n_cols_plot)

width_page = 12
height_page = 1.5 * n_rows_plot #Inches
format = file_ext(output_file_name)

if(format == 'pdf' | format == 'PDF'){
  pdf(output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
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


