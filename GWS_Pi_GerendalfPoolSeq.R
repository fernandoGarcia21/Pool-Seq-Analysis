setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
library(RColorBrewer)

generate_intercalated_colors <- function(n, color1 = "#80cdc1", color2 = "#dfc27d") {
  # Create two sequences of colors, each with a length sufficient to cover n
  seq1 <- rep(color1, ceiling(n / 2))
  seq2 <- rep(color2, floor(n / 2))
  
  # Interleave the two sequences
  interleaved_colors <- character(n)
  interleaved_colors[seq(1, n, by = 2)] <- seq1[1:ceiling(n / 2)]
  interleaved_colors[seq(2, n, by = 2)] <- seq2[1:floor(n / 2)]
  
  return(interleaved_colors)
}


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
    cat("\nNo files found in the specified folder.\n")
    merged_df = NULL # Or however you want to handle the case with no files
  }
  return(merged_df)
}



#Directory where the diversity statistics are found
#files_dir_diversity_path = 'Grenedalf/TajimasD_100K';
files_dir_diversity_outliers_path = 'Grenedalf/FST_PI_20K';

files_dir_diversity_path = 'Grenedalf/WindowDiversity_20K';

#output_file_name = 'Europe_PiRatioFST_Grenedalf.pdf';
output_file_name = 'Pi_Grenedalf.NonBiased.jpg';

focus_population = 'RED'
area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['SM-BOO-low-merged']] = 'NA'
area_focus_population[['PLY-high']] = 'NA'
area_focus_population[['VEN']] = 'EU'
area_focus_population[['EAS-low']] = 'EU'

output_file_name = paste(focus_population, output_file_name, sep = '.')

# List all files in the specified folder with the .diversity.csv extension
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)
merged_diversity_df = read_and_merge_files(file_list_diversity, c(focus_population), '\t')


#Read outliers file OAK.SM-BOO-low-merged.WindowOutliers
file_list_diversity_outliers = list.files(path = files_dir_diversity_outliers_path, pattern = paste0("^",focus_population,".*\\.WindowOutliers\\.NonBiased\\.txt$"), full.names = TRUE)
#file_list_diversity_outliers = list.files(path = files_dir_diversity_outliers_path, pattern = paste0("^",focus_population,".*\\.sharedOutlierWindows\\.txt$"), full.names = TRUE)
if(length(file_list_diversity_outliers) > 0){
  diversity_outliers_df = read.table(file_list_diversity_outliers[1], header = TRUE, sep = '\t', check.names = FALSE)
}else{
  diversity_outliers_df = data.frame()
}

print('Finish reading datasets')
print(paste('Windows in the Diversity dataset:', dim(merged_diversity_df)[1]))


#Subset the dataframe to get only the pi of the focus population
tmp_col_pi = paste( focus_population, 'theta_pi' , sep = '.')
tmp_stats_pop_df = merged_diversity_df[,c('chrom','start','end',tmp_col_pi)]
colnames(tmp_stats_pop_df) =  c('chrom','start','end','pi')

#Remove NA rows
tmp_stats_pop_df = na.omit(tmp_stats_pop_df) #Remove rows with NA values 

#Remove the whole diversity DF to release memmory
rm(merged_diversity_df)



#Prepare the metadata to plot the X axis

#Extract the LG number from the LG name
tmp_stats_pop_df$LG = sapply(strsplit(tmp_stats_pop_df$chrom, "_"), function(x) x[1])
tmp_stats_pop_df$LGN = as.numeric(substr(tmp_stats_pop_df$LG,3,length(tmp_stats_pop_df$LG)))


#Sort the matrix by the LG number and window number
tmp_stats_pop_df = tmp_stats_pop_df[order(tmp_stats_pop_df$LGN, tmp_stats_pop_df$start),]


#Compute the cumulative of the windows to plot each value after the other
pi_df_sorted_final = tmp_stats_pop_df %>% 
  # Compute chromosome size
  group_by(LGN) %>% 
  summarise(LG_len=max(start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(LG_len)-LG_len) %>%
  select(-LG_len) %>%
  
  # Add this info to the initial dataset
  left_join(tmp_stats_pop_df, ., by=c("LGN"="LGN")) %>%
  
  # Add a cumulative position of each SNP
  arrange(LGN, start) %>%
  mutate( BPcum=start+tot)


# Determine the x position of the axis labels
axis_pi_df <- pi_df_sorted_final %>% 
  group_by(LGN) %>% 
  summarize(start=min(BPcum),
            end=max(BPcum))
#Add the gray white intercalated
axis_pi_df$color <- rep(c(1,0),ceiling(nrow(axis_pi_df)/2))[1:nrow(axis_pi_df)]

#Estimate the center of the chromosome
ref_end_previouslg = append(0, axis_pi_df$end+0.0001)
ref_end_previouslg = ref_end_previouslg[1:length(ref_end_previouslg)-1]
axis_pi_df$start = ref_end_previouslg
axis_pi_df$center = axis_pi_df$start + ((axis_pi_df$end - axis_pi_df$start) / 2)

#Generate a list of intercalated colours to distinguish the chromosomes
n = length(unique(pi_df_sorted_final$LGN))
color_array = generate_intercalated_colors(n, color1 = "#80cdc1", color2 = "#dfc27d")
unique_lgn = unique(pi_df_sorted_final$LGN)
listColourCode = setNames(as.list(color_array), unique_lgn) 

# Create a new 'color' column by looking up the color for each LGN value
pi_df_sorted_final$color = listColourCode[as.character(pi_df_sorted_final$LGN)]



#Estimate the average Pi Every 100 windows

pi_df_averaged <- pi_df_sorted_final %>%
  group_by(LGN) %>%
  mutate(row_number_in_group = row_number()) %>%
  group_by(LGN, group = ceiling(row_number_in_group / 100)) %>%
  summarise(average_pi = mean(pi, na.rm = TRUE),
            min_BPcum = min(BPcum, na.rm = TRUE),
            max_BPcum = max(BPcum, na.rm = TRUE),
            n_rows = n(),
            .groups = 'drop') %>%
  ungroup() %>%
  select(LGN, group, average_pi, min_BPcum, max_BPcum, n_rows)


color_array = generate_intercalated_colors(n, color1 = "#018571", color2 = "#a6611a")
listColourCode = setNames(as.list(color_array), unique_lgn) 

# Create a new 'color' column by looking up the color for each LGN value
pi_df_averaged$color = listColourCode[as.character(pi_df_averaged$LGN)]

#Add the coordinate of the X and Y where the segment ends
pi_df_averaged <- pi_df_averaged %>%
  group_by(LGN) %>% # Optional: If you want to calculate Y_end within groups
  mutate(X_end = lead(min_BPcum, n = 1, default = NA)) %>%
  mutate(Y_end = lead(average_pi, n = 1, default = NA)) %>%
  ungroup()



#Update the coordinates of the outliers DF to the cumulative coordinates on the X axis for plotting

coordinates_outliers_df <- pi_df_sorted_final %>%
  inner_join(diversity_outliers_df, by = c("chrom" = "chrom", "start" = "start"))
coordinates_outliers_df = coordinates_outliers_df[,c('chrom','start','BPcum','LGN')]

#coordinates_outliers_df$X_end = coordinates_outliers_df$BPcum + window_size


#Change colour to the Pi windows that are part of the outliers
pi_df_sorted_final$group = 'all'
#pi_df_sorted_final[pi_df_sorted_final$BPcum %in% coordinates_outliers_df$BPcum, ]$color = '#f30057'
#pi_df_sorted_final[pi_df_sorted_final$BPcum %in% coordinates_outliers_df$BPcum, ]$group = 'outlier'


# Calculate the 95th and 5th percentile of the pi distribution
tmp_percentile = 99
percentile_upper = quantile(pi_df_sorted_final$pi, probs = tmp_percentile/100, na.rm = TRUE)
percentile_lower = quantile(pi_df_sorted_final$pi, probs = 1 - tmp_percentile/100, na.rm = TRUE)

pi_df_sorted_final[pi_df_sorted_final$pi < percentile_lower, ]$color = '#f30057'
pi_df_sorted_final[pi_df_sorted_final$pi < percentile_lower, ]$group = 'outlier'


pltPi = ggplot()+
  #Draw veritical bars representing the outlier windows
  #geom_rect(data=coordinates_outliers_df, mapping=aes(xmin=BPcum, xmax=X_end, ymin=-4, ymax=1), fill='#f30057', show.legend = FALSE)+
  
  # Show all points
  geom_point(data= pi_df_sorted_final[pi_df_sorted_final$group == 'all',], aes(x=BPcum, y=pi, color = color), size = 0.1, alpha=0.5) +
  # Show all points
  geom_point(data= pi_df_sorted_final[pi_df_sorted_final$group == 'outlier',], aes(x=BPcum, y=pi, color = color), size = 0.5, alpha=0.5) +
  #Line that represents the average statistic for a given set of windows
  geom_segment(data = pi_df_averaged, aes(x=min_BPcum, y=average_pi, xend = X_end, yend = Y_end, color = color))+
  
  # Add the horizontal red line at the 95th percentile
  #geom_hline(yintercept = percentile_lower, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
  # custom X axis:
  scale_x_continuous( label = axis_pi_df$LGN, breaks= axis_pi_df$center , expand = c(0.01,0)) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
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
    plot.title = element_text(hjust = 0.5, face='bold', size=12,
                              margin = margin(t = 0, r = 0, b = -50, l = 0)),
    axis.title=element_text(size=10,face="bold"),
    axis.text.y=element_text(size=8),
    axis.text.x=element_text(size=12),
    plot.margin = margin( 0,0.0,0.0,0.0, "cm")
    
  )+
  ylim(0, 0.06)+ #For Oak and Red
  #ylim(0, 0.02)+ #for Venice
  labs(y = "Theta Pi", x="") +
  ggtitle(paste(focus_population, 'Theta Pi'))+
  #geom_rect(data=axis_pi_df, mapping=aes(xmin=start, xmax=end, ymin=-4, ymax=1, fill=as.factor(color)), alpha=0.15, show.legend = FALSE)+
  #scale_fill_manual(values = c('1'="gray", '0'="white"), name=NULL)+
  guides(fill = "none")




width_page = 15
height_page = 4 
format = file_ext(output_file_name)

if(format == 'pdf' | format == 'PDF'){
  pdf(output_file_name, width = width_page, height = height_page)
  print(pltPi)
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    ggsave(output_file_name, 
           plot = pltPi, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')



