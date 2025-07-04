setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
library(RColorBrewer)
library(ggvenn)
library(VennDiagram)

read_and_merge_files = function(p_file_list, p_sep){
  
  # Initialize an empty list to store the data frames
  data_frames_list <- list()
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



##############################################################################################
# Plots genome coordinates and chromosomes along the x axis vs a statistic like FST on the Y axis
##############################################################################################
generate_manhattan_plot = function (p_statistic_df, p_focus_population, p_pop_vs, p_suffix_stat, p_stat_name, p_percentile, p_outliers_df, is_above_percentile, p_max_y){
  #Subset the dataframe to get only the p_suffix_stat statistics for the current vs population
  tmp_col_stat = ''
  if(is.null(p_pop_vs)){
    tmp_col_stat = paste0( p_focus_population, p_suffix_stat)
    
  }else{
    tmp_col_stat = paste0( p_focus_population, ':', p_pop_vs, p_suffix_stat)
  }
  
  
  tmp_stats_pop_df = p_statistic_df[,c('chrom','start','end',tmp_col_stat)]
  tmp_stats_pop_df = na.omit(tmp_stats_pop_df) #Remove rows with NA values 
  colnames(tmp_stats_pop_df) =  c('chrom','start','end', 'stat')
  
  
  #Prepare the metadata to plot the X axis
  
  #Extract the LG number from the LG name
  tmp_stats_pop_df$LG = sapply(strsplit(tmp_stats_pop_df$chrom, "_"), function(x) x[1])
  tmp_stats_pop_df$LGN = as.numeric(substr(tmp_stats_pop_df$LG,3,length(tmp_stats_pop_df$LG)))
  
  
  #Sort the matrix by the LG number and window number
  tmp_stats_pop_df = tmp_stats_pop_df[order(tmp_stats_pop_df$LGN, tmp_stats_pop_df$start),]
  
  
  #Compute the cumulative of the windows to plot each value after the other
  stat_df_sorted_final = tmp_stats_pop_df %>% 
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
  axis_stat_df <- stat_df_sorted_final %>% 
    group_by(LGN) %>% 
    summarize(start=min(BPcum),
              end=max(BPcum))
  #Add the gray white intercalated
  axis_stat_df$color <- rep(c(1,0),ceiling(nrow(axis_stat_df)/2))[1:nrow(axis_stat_df)]
  
  #Estimate the center of the chromosome
  ref_end_previouslg = append(0, axis_stat_df$end+0.0001)
  ref_end_previouslg = ref_end_previouslg[1:length(ref_end_previouslg)-1]
  axis_stat_df$start = ref_end_previouslg
  axis_stat_df$center = axis_stat_df$start + ((axis_stat_df$end - axis_stat_df$start) / 2)
  
  #Generate a list of intercalated colours to distinguish the chromosomes
  n = length(unique(stat_df_sorted_final$LGN))
  color_array = generate_intercalated_colors(n, color1 = "#80cdc1", color2 = "#dfc27d")
  unique_lgn = unique(stat_df_sorted_final$LGN)
  listColourCode = setNames(as.list(color_array), unique_lgn) 
  
  # Create a new 'color' column by looking up the color for each LGN value
  stat_df_sorted_final$color = listColourCode[as.character(stat_df_sorted_final$LGN)]
  
  
  
  #Estimate the average stat Every 100 windows
  
  stat_df_averaged <- stat_df_sorted_final %>%
    group_by(LGN) %>%
    mutate(row_number_in_group = row_number()) %>%
    group_by(LGN, group = ceiling(row_number_in_group / 100)) %>%
    summarise(average_stat = mean(stat, na.rm = TRUE),
              min_BPcum = min(BPcum, na.rm = TRUE),
              max_BPcum = max(BPcum, na.rm = TRUE),
              n_rows = n(),
              .groups = 'drop') %>%
    ungroup() %>%
    select(LGN, group, average_stat, min_BPcum, max_BPcum, n_rows)
  
  
  color_array = generate_intercalated_colors(n, color1 = "#018571", color2 = "#a6611a")
  listColourCode = setNames(as.list(color_array), unique_lgn) 
  
  # Create a new 'color' column by looking up the color for each LGN value
  stat_df_averaged$color = listColourCode[as.character(stat_df_averaged$LGN)]
  
  #Add the coordinate of the X and Y where the segment ends
  stat_df_averaged <- stat_df_averaged %>%
    group_by(LGN) %>% # Optional: If you want to calculate Y_end within groups
    mutate(X_end = lead(min_BPcum, n = 1, default = NA)) %>%
    mutate(Y_end = lead(average_stat, n = 1, default = NA)) %>%
    ungroup()
  
  
  
  # Calculate the 95th and 5th percentile of the stat distribution
  tmp_percentile = p_percentile
  percentile_upper = quantile(stat_df_sorted_final$stat, probs = tmp_percentile/100, na.rm = TRUE)
  percentile_lower = quantile(stat_df_sorted_final$stat, probs = 1 - tmp_percentile/100, na.rm = TRUE)
  
  #Colour the outliers based on the 9th percentile of the  distribution
  if(is_above_percentile){
    stat_df_sorted_final[stat_df_sorted_final$stat > percentile_upper, ]$color = '#fcb603'
    stat_df_sorted_final[stat_df_sorted_final$stat > percentile_upper, ]$group = 'outlier'
  }else{
    stat_df_sorted_final[stat_df_sorted_final$stat < percentile_lower, ]$color = '#fcb603'
    stat_df_sorted_final[stat_df_sorted_final$stat < percentile_lower, ]$group = 'outlier'
  }
  
  
  #Update the coordinates of the outliers DF to the cumulative coordinates on the X axis for plotting
  coordinates_outliers_df <- stat_df_sorted_final %>%
    inner_join(p_outliers_df, by = c("chrom" = "chrom", "start" = "start"))
  coordinates_outliers_df = coordinates_outliers_df[,c('chrom','start','BPcum','LGN')]
  
  
  #Change colour to the windows that are part of the outliers dataframe
  stat_df_sorted_final$group = 'all'
  #Colour the outliers based on the fst vs pi windows
  stat_df_sorted_final[stat_df_sorted_final$BPcum %in% coordinates_outliers_df$BPcum, ]$color = '#f30057'
  stat_df_sorted_final[stat_df_sorted_final$BPcum %in% coordinates_outliers_df$BPcum, ]$group = 'candidate'
  
  
  pltStat = ggplot()+
    #Draw veritical bars representing the outlier windows
    #geom_rect(data=coordinates_outliers_df, mapping=aes(xmin=BPcum, xmax=X_end, ymin=-4, ymax=1), fill='#f30057', show.legend = FALSE)+
    
    # Show all points
    geom_point(data= stat_df_sorted_final[stat_df_sorted_final$group == 'all',], aes(x=BPcum, y=stat, color = color), size = 0.1, alpha=0.5) +
    # Show outlier points
    geom_point(data= stat_df_sorted_final[stat_df_sorted_final$group == 'outlier',], aes(x=BPcum, y=stat, color = color), size = 0.5, alpha=0.5) +
    # Show candidate points, outliers with high fst and low diversity
    geom_point(data= stat_df_sorted_final[stat_df_sorted_final$group == 'candidate',], aes(x=BPcum, y=stat, color = color), size = 1.5, alpha=0.5) +
    #Line that represents the average statistic for a given set of windows
    geom_segment(data = stat_df_averaged, aes(x=min_BPcum, y=average_stat, xend = X_end, yend = Y_end, color = color))+
    
    # Add the horizontal red line at the 95th percentile
    #geom_hline(yintercept = percentile_lower, color = "darkgray", linetype = "dashed", linewidth = 0.8) +
    # custom X axis:
    scale_x_continuous( label = axis_stat_df$LGN, breaks= axis_stat_df$center , expand = c(0.01,0)) +
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
                                margin = margin(t = 0, r = 0, b = -3, l = 0)),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=8),
      axis.text.x=element_text(size=12),
      plot.margin = margin( 0.3,0.0,0.0,0.0, "cm")
      
    )+
    ylim(-0.01, p_max_y)+
    labs(y = p_stat_name, x="") +
    ggtitle(paste(p_stat_name, p_focus_population, 'vs', p_pop_vs))+
    #geom_rect(data=axis_stat_df, mapping=aes(xmin=start, xmax=end, ymin=-4, ymax=1, fill=as.factor(color)), alpha=0.15, show.legend = FALSE)+
    #scale_fill_manual(values = c('1'="gray", '0'="white"), name=NULL)+
    guides(fill = "none")
  
  return (pltStat)
}


keyword = 'LowPiNoThreshold'
#Directory where the FST windows are found, we just need this to know how many windows were used
files_dir_fst_path = 'Grenedalf/WindowFST_10K';
files_dir_diversity_path = 'Grenedalf/WindowDiversity_10K';

#output_file_name = 'Europe_PiRatioFST_Grenedalf.pdf';
output_file_name = paste0(keyword,'Shared_FSTLowPiRegions_Grenedalf10k.NonBiased.jpg');

fst_percentile = 95 #FST Threshold for North America is 95, and 99 for Venice?
pi_percentile = 95
#focus_population = 'Maine'
focus_population = 'VEN'

list_pops_vs_area = list()
list_pops_vs_area[['OAK']] = c('SM-BOO-high-merged', 'SM-BOO-low-merged', 'PBY_BLMP1', 'PTY_BLMP2', 'PLY-high')
list_pops_vs_area[['RED']] = c('SM-BOO-high-merged', 'SM-BOO-low-merged', 'PBY_BLMP1', 'PTY_BLMP2', 'PLY-high')
list_pops_vs_area[['PLY-high']] = c('PLY-high')
list_pops_vs_area[['SM-BOO-high-merged']] = c('SM-BOO-high-merged')
list_pops_vs_area[['SM-BOO-low-merged']] = c('SM-BOO-low-merged')
list_pops_vs_area[['PTY_BLMP2']] = c('PTY_BLMP2')
list_pops_vs_area[['VEN']] = c('EAS-low', 'T_NC_NO-high', 'Fr_Crab_Low', 'EAS-high', '3GAL')

area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['PLY-high']] = 'NA'
area_focus_population[['Maine']] = 'NA'
area_focus_population[['SM-BOO-low-merged']] = 'NA'
area_focus_population[['VEN']] = 'EU'
area_focus_population[['EAS-low']] = 'EU'

output_file_name = paste(focus_population, output_file_name, sep = '.')


# List all files in the specified folder with the .diversity.csv extension
file_list_fst = list.files(path = files_dir_fst_path, pattern = paste0("\\.",focus_population,".*\\.fst\\.csv$"), full.names = TRUE)
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)

merged_fst_df = read_and_merge_files(file_list_fst, '\t')
merged_diversity_df = read_and_merge_files(file_list_diversity, '\t')

# Merge the dataframes by the first three columns
merged_df_all_windows = merge(merged_fst_df, merged_diversity_df, by = c("chrom", "start", "end"))

# Create a unique identifier for all windows DF
merged_df_all_windows = merged_df_all_windows %>%
  mutate(unique_id = paste(chrom, start, end, sep = "_"))


list_pops_vs = list_pops_vs_area[[focus_population]]
stats_pop_vs_outliers_df = data.frame()
tmp_list_pi_thresholds = list()
for(pop_vs in list_pops_vs){
  
  #Subset the dataframe to get only the fst and pi statistics for the current vs population
  tmp_col_fst = paste0( focus_population, ':', pop_vs, '.fst')
  tmp_cols_pi = paste0( c(focus_population, pop_vs), '.theta_pi')
  
  tmp_stats_pop_vs = merged_df_all_windows[,c('chrom','start','end', 'unique_id',tmp_col_fst,tmp_cols_pi)]
  tmp_stats_pop_vs = na.omit(tmp_stats_pop_vs) #Remove rows with NA values 
  colnames(tmp_stats_pop_vs) =  c('chrom','start','end', 'unique_id', 'fst','pi_focus', 'pi_vs')
  
  
  #Identify the FST outliers based on the FST percentil
  tmp_fst_percentile = fst_percentile
  percentile_fst = quantile(tmp_stats_pop_vs$fst, probs = tmp_fst_percentile/100, na.rm = TRUE)
  
  #Identify the outliers in the candidate source population based on the pi percentile
  percentile_pi_vs = quantile(tmp_stats_pop_vs$pi_vs, probs = 1-(pi_percentile/100), na.rm = TRUE)
  
  #Identify the outliers in the introduced population based on the pi percentile (e.g. top 5% lowest)
  percentile_pi_focus = quantile(tmp_stats_pop_vs$pi_focus, probs = 1-(pi_percentile/100), na.rm = TRUE)
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$ For extra validation $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  tmp_list_pi_thresholds[[focus_population]] = percentile_pi_focus
  tmp_list_pi_thresholds[[pop_vs]] = percentile_pi_vs
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  #Identify FST outliers
  tmp_stats_pop_vs_outliers = tmp_stats_pop_vs[tmp_stats_pop_vs$fst >= percentile_fst, ]
  
  print(paste('^^^^^^^^^^^^^^^^^^^^^^^^ Percentage of outliers: ', dim(tmp_stats_pop_vs_outliers)[1] / dim(tmp_stats_pop_vs) [1]))
  
  #If FST outliers were found, identify the fst outliers where the nucleotide diversity is reduced in the focus population
  if(dim(tmp_stats_pop_vs_outliers)[1] > 0){
    tmp_stats_pop_vs_outliers$type = 'FST_Only'
    tmp_stats_pop_vs_outliers$group = pop_vs

    tmp_stats_pop_vs_outliers[tmp_stats_pop_vs_outliers$pi_focus < tmp_stats_pop_vs_outliers$pi_vs
                     , ]$type = 'FST_and_Pi'
    
    #Append the temporal outliers dataframe to final dataframe
    stats_pop_vs_outliers_df = rbind(stats_pop_vs_outliers_df, tmp_stats_pop_vs_outliers)
    
  }else{
    print('No outliers found')
  }
  
}#For iterating over vs populations


########################################################################
# Print the shared outliers
#########################################################################
if(dim(stats_pop_vs_outliers_df)[1] > 0){
  
  #Write the total outliers for all populations
  tmp_outlier_file_name = paste(focus_population, keyword, 'AllFstOutlierLowPiRegionsWindows.txt', sep = '.')
  write.table(stats_pop_vs_outliers_df, file = tmp_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  
  #Write the true outliers for all populations; only fst_and_pi outliers
  tmp_outlier_file_name = paste(focus_population, keyword, 'TrueOutlierLowPiRegionsWindows.txt', sep = '.')
  write.table(stats_pop_vs_outliers_df[stats_pop_vs_outliers_df$type == 'FST_and_Pi', ], file = tmp_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  
   #########################################################################
  # PLOT VENN DIAGRAM 
  ##########################################################################
  
  #Identify only unique windows (removing, fst, pi, group and type information)
  diversity_outliers_df = stats_pop_vs_outliers_df[stats_pop_vs_outliers_df$type == 'FST_and_Pi', ]
  # Create a unique identifier for all windows DF
  #diversity_outliers_df = diversity_outliers_df %>%
  #  mutate(unique_id = paste(chrom, start, end, sep = "_"))
  
  
  
  # Function to create the list for ggvenn
  create_venn_list = function(group_name) {
    diversity_outliers_df %>%
      filter(group == group_name) %>%
      pull(unique_id)
  }
  
  # Get unique groups
  groups = unique(diversity_outliers_df$group)
  
  # Create the list
  venn_list <- lapply(groups, create_venn_list)
  names(venn_list) <- groups # Name the list elements
  
  # Print the list (for debugging)
  #print(venn_list)
  
  #Identify the elements that are common to the list
  shared_outliers = Reduce(intersect, venn_list)
  
  
  
  myCol = brewer.pal(length(venn_list), "Set2")
  
  # Chart
  venn.diagram(
    x = venn_list,
    filename = output_file_name,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 1080 , 
    width = 1080 , 
    resolution = 300,
    compression = "lzw",
    main = paste(focus_population, 'Outlier windows'),
    main.cex = 0.6,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = myCol,
    #col = 'darkgray',
    
    # Numbers
    cex = .6,
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
    
    #margin to give space for the names
    #margin = 0.3
  )
  
  
  
  print('Output file created')
  
  
}



########################################################################
# Print the shared outliers and estimate the probability of being shared
#########################################################################
if(length(shared_outliers) > 0){
  
  shared_diversity_outliers_df = diversity_outliers_df[diversity_outliers_df$unique_id %in% shared_outliers, c('chrom','start','end')]
  shared_diversity_outliers_df = unique(shared_diversity_outliers_df)
  #Export the shared outlier windows
  tmp_shared_outlier_file_name = paste(focus_population, keyword, 'sharedFstLowPiRegionsOutlierWindows.txt', sep = '.')
  write.table(shared_diversity_outliers_df, file = tmp_shared_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  

  #Identify the universe of outliers among all populations
  outliers_universe_df = diversity_outliers_df[, c('chrom','start','end')]
  outliers_universe_df = unique(outliers_universe_df)
  
  #Count how many outliers were identified for each versus population
  outliers_by_population = diversity_outliers_df %>%
    group_by(group) %>%
    summarise(num_outliers = n())
  
  #Percentage of shared outliers of of the total
  print(paste("The proporcion of shared outliers out of the total unique: ", length(shared_outliers)/nrow(outliers_universe_df)))
  
  #Universe of outliers
  print(paste("The total number number of unique outlier windows is: ", nrow(outliers_universe_df)))
  
  #expected number of Shared outliers by chance
  shared_expectation = prod(outliers_by_population$num_outliers) / (nrow(outliers_universe_df)^(nrow(outliers_by_population)-1))
  print(paste("Number of outliers expected shared by chance: ", shared_expectation))
  
  #Estimate the p-value
  # Merge the dataframes by the first three columns
  # Create a unique identifier for all windows DF
  merged_outliers_universe_df <- outliers_universe_df %>%
    mutate(unique_id = paste(chrom, start, end, sep = "_"))
  
  
  # Add new columns with FALSE values and add TRUE to the windows that are true in each group of outliers
  for (group in outliers_by_population$group) {
    merged_outliers_universe_df[[group]] = FALSE
    tmp_outliers_group = diversity_outliers_df[diversity_outliers_df$group == group,]
    merged_outliers_universe_df[merged_outliers_universe_df$unique_id %in% tmp_outliers_group$unique_id, group] = TRUE
  }
  
  #Initialize parameters for permutation
  observed_shared_outliers = length(shared_outliers) # The number of shared outliers you observed
  number_groups = nrow(outliers_by_population)
  group_cols = outliers_by_population$group
  count_chance_shared_occurrences = 0
  perm_replicates = 10000
  
  
  for(i in 1:perm_replicates){
    tmp_shuffled_windows = merged_outliers_universe_df
    
    # Shuffle the TRUE and FALSE labels within each group column independently
    for(group in group_cols){
      tmp_shuffled_windows[[group]] = sample(tmp_shuffled_windows[[group]], replace = FALSE)
    }
    
    # Identify windows that are TRUE in ALL group columns in this shuffled dataset
    shared_outliers_permutation = tmp_shuffled_windows %>%
      filter(if_all(all_of(group_cols), ~ . == TRUE))
    
    num_shared_permutation = nrow(shared_outliers_permutation)
    
    
    # Check if the number of shared outliers in this permutation is >= the observed number
    if(num_shared_permutation >= observed_shared_outliers){
      print(paste('Shared in permutation', i, ":", num_shared_permutation) )
      count_chance_shared_occurrences = count_chance_shared_occurrences + 1
    }
  }
  
  # Calculate the p-value
  p_value = count_chance_shared_occurrences / perm_replicates
  print(paste("Estimated p-value:", p_value))
  
  
  #Resolution of p-value: The smallest p-value you can theoretically resolve with N permutations is 1/(N+1).
  #With 1,000 permutations, the smallest p-value you can resolve (other than 0) is approximately 0.001.
  #With 10,000 permutations, the smallest p-value you can resolve is approximately 0.0001.
}




#################################################################################
# Plot FST Manhattan with the windows and outliers 
#################################################################################

#For loop iterating over the populations in the list_pops_vs
list_fst_plt = list()
tmp_colnames_df = colnames(merged_fst_df)
list_pops_vs = list_pops_vs_area[[focus_population]]
aux_above_percentile = TRUE
for(pop_vs in list_pops_vs){
  
  tmp_plt_fst = generate_manhattan_plot(merged_fst_df, focus_population, pop_vs, '.fst', 'FST', fst_percentile, shared_diversity_outliers_df, aux_above_percentile, 1)
  
  list_fst_plt[[paste0('FST_',pop_vs)]] = tmp_plt_fst
}


#Plot a grid the FST Manhattans

n_cols_plot = 1
n_rows_plot = length(list_fst_plt)

width_page = 12
height_page = 3.5 * n_rows_plot #Inches
format = file_ext(output_file_name)

tmp_output_file_name = paste('Manhattan', output_file_name, sep = '.')

if(format == 'pdf' | format == 'PDF'){
  pdf(tmp_output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_fst_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_fst_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
    ggsave(tmp_output_file_name, 
           plot = combined_page, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')





#################################################################################
# Plot Diversity Manhattan with the windows and outliers 
#################################################################################

#For loop iterating over the populations in the list_pops_vs
list_pi_plt = list()
tmp_colnames_df = colnames(merged_diversity_df)
list_pops_vs = list_pops_vs_area[[focus_population]]
list_pops_vs = append(focus_population, list_pops_vs)
aux_above_percentile = FALSE
for(pop_vs in list_pops_vs){
  
  tmp_plt_pi = generate_manhattan_plot(merged_diversity_df, pop_vs, NULL, '.theta_pi', 'Theta pi', pi_percentile, shared_diversity_outliers_df, aux_above_percentile, 0.06) #Max pi for NA must be 0.15 
  
  list_pi_plt[[paste0('Pi_',pop_vs)]] = tmp_plt_pi
}


#Plot a grid the FST Manhattans

n_cols_plot = 1
n_rows_plot = length(list_pi_plt)

width_page = 12
height_page = 3.5 * n_rows_plot #Inches
format = file_ext(output_file_name)

tmp_output_file_name = paste('ManhattanPi', output_file_name, sep = '.')

if(format == 'pdf' | format == 'PDF'){
  pdf(tmp_output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_pi_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_pi_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
    ggsave(tmp_output_file_name, 
           plot = combined_page, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')



#################################################################################
# Distribution of distances between consecutive windows
#################################################################################

shared_diversity_outliers_df = diversity_outliers_df[diversity_outliers_df$unique_id %in% shared_outliers, c('chrom','start','end')]
shared_diversity_outliers_df = unique(shared_diversity_outliers_df)

#Count how many outliers were identified for each versus population
outliers_by_population = diversity_outliers_df %>%
  group_by(group) %>%
  summarise(num_outliers = n())

group_cols = outliers_by_population$group


cat("Original data frame (first few rows):\n")
print(head(shared_diversity_outliers_df))
cat("\n")


# Calculate distances
distances_df <- shared_diversity_outliers_df %>%
  dplyr::group_by(chrom) %>%
  dplyr::arrange(start) %>%
  dplyr::mutate(
    next_start = dplyr::lead(start),
    distance_bp = next_start - end
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(distance_bp))


# --- Plotting the histogram using ggplot2 ---
distances_df$log_distance_bp = log(distances_df$distance_bp, base = 10)

plt_actual_distances = ggplot(distances_df, aes(x = log_distance_bp))+
  geom_histogram(binwidth = 0.5, fill = "#80cdc1", color = "black", linewidth = 0.2, alpha = 0.5,
                 center = 0.25) + # <--- Add center = 0.25 here
  labs(
    title = paste(focus_population, "- Distances Between Consecutive Genomic Windows"),
    x = "Distance log10(bp)",
    y = "Frequency"
  ) +
  theme_minimal() +
  scale_x_continuous(
    limits = c(0, 10),                 # Set the exact range of the x-axis
    breaks = seq(0, 10, by = 1)        # Set tick marks at every integer from 0 to 10
  ) +
  theme(
    legend.position="none",
    legend.title = element_text(size = 8, face = "bold"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(size = 0.1),
    plot.title = element_text(hjust = 0.5, face='bold', size=12,
                              margin = margin(t = 0, r = 0, b = -3, l = 0)),
    axis.title=element_text(size=10,face="bold"),
    axis.text.y=element_text(size=8),
    axis.text.x=element_text(size=12),
    plot.margin = margin( 0.3,0.0,0.0,0.0, "cm")
  )



############## Permutation test of distances 

# --- 2. Calculate the Observed Test Statistic ---
observed_median_distance <- median(distances_df$log_distance_bp)
cat("Observed Median Distance:", observed_median_distance, "\n")

perm_replicates = 10000
distance_threshold = log(10000, base = 10)

proportion_observed_below_threshold = length(distances_df$log_distance_bp[distances_df$log_distance_bp <= distance_threshold]) / nrow(distances_df)
cat("Observed distance below", distance_threshold, ": ", proportion_observed_below_threshold, "\n")


#Identify the universe of outliers among all populations
outliers_universe_df = diversity_outliers_df[, c('chrom','start','end')]
outliers_universe_df = unique(outliers_universe_df)

#Estimate the p-value
# Merge the dataframes by the first three columns
# Create a unique identifier for all windows DF
merged_outliers_universe_df <- outliers_universe_df %>%
  mutate(unique_id = paste(chrom, start, end, sep = "_"))


# Add new columns with FALSE values and add TRUE to the windows that are true in each group of outliers
for (group in outliers_by_population$group) {
  merged_outliers_universe_df[[group]] = FALSE
  tmp_outliers_group = diversity_outliers_df[diversity_outliers_df$group == group,]
  merged_outliers_universe_df[merged_outliers_universe_df$unique_id %in% tmp_outliers_group$unique_id, group] = TRUE
}


perm_distances_list = list()
perm_prop_below_threshold = c()
for(i in 1:perm_replicates){
  tmp_shuffled_windows = merged_outliers_universe_df
  
  # Shuffle the TRUE and FALSE labels within each group column independently
  for(group in group_cols){
    tmp_shuffled_windows[[group]] = sample(tmp_shuffled_windows[[group]], replace = FALSE)
  }
  
  # Identify windows that are TRUE in ALL group columns in this shuffled dataset
  shared_outliers_permutation = tmp_shuffled_windows %>%
    filter(if_all(all_of(group_cols), ~ . == TRUE))
  
  # Calculate distances
  distances_permutation_df <- shared_outliers_permutation %>%
    dplyr::group_by(chrom) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(
      next_start = dplyr::lead(start),
      distance_bp = next_start - end
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(distance_bp))
  
  
  #Append the estimates of distance to the list
  distances_permutation_df$log_distance_bp = log(distances_permutation_df$distance_bp, base = 10)
  perm_distances_list[[i]] = distances_permutation_df$log_distance_bp
  
  proportion_perm_below_threshold = length(distances_permutation_df$log_distance_bp[distances_permutation_df$log_distance_bp <= distance_threshold]) / nrow(distances_permutation_df)
  perm_prop_below_threshold = append(perm_prop_below_threshold, proportion_perm_below_threshold)
}

print(perm_prop_below_threshold)

# This p-value estimate is based on the number of distances below the threshold. To reject the null hypothesis that shuffled
# labels will produce distances smaller or equal to the actual data
p_value_less_threshold <- (sum(proportion_observed_below_threshold <= perm_prop_below_threshold ) + 1) / (perm_replicates + 1)

cat("P-value (testing if proportion of observed distances below ", distance_threshold, " is significantly smaller):", p_value_less_threshold, "\n")
cat("\n")


# --- 3. Calculate Test Statistics for Each Permutation ---
permuted_median_distances = sapply(perm_distances_list, median)

# Add 1 to numerator and denominator to handle cases where observed_median_distance is smaller than all permuted_median_distances
# This provides a more conservative and appropriate p-value for permutation tests when R=0.
p_value <- (sum(permuted_median_distances <= observed_median_distance) + 1) / (perm_replicates + 1)

cat("P-value (testing if observed median is significantly smaller):", p_value, "\n")
cat("\n")


# using sprintf() for consistent formatting, even for whole numbers like 7 -> "7.0"
formatted_median_distance = sprintf("%.1f", observed_median_distance)


# Option B: Histogram of Permuted Statistics with Observed Statistic Line
# This plot shows where your observed median falls relative to the null distribution of medians.
plt_median_distances = ggplot(data.frame(Permuted_Medians = permuted_median_distances), aes(x = Permuted_Medians)) +
  geom_histogram(binwidth = 0.1, fill = "#a6611a", color = "black", linewidth = 0.2, alpha = 0.5) + # Adjust binwidth as needed
  geom_vline(aes(xintercept = observed_median_distance), color = "red", linetype = "dashed", linewidth = 1) +
  # Add the text annotation
  annotate(
    "text",
    x = observed_median_distance + 0.5,  # X-position is the line itself
    y = Inf,                       # Place it at the top of the plot (Inf for infinity)
    label = formatted_median_distance, # Use the formatted value
    color = "red",
    vjust = 1.5,                   # Adjust vertical position down from Inf (1.5 is a common good starting point)
    hjust = -0.1,                  # Adjust horizontal position (e.g., -0.1 to nudge right)
    size = 3.5,                    # Adjust text size as needed
    fontface = "bold"              # Make text bold if desired
  ) +
  labs(
    title = paste(focus_population, "- Distribution of Permuted Median Distances"),
    x = "Median Distance log10(bp)",
    y = "Frequency"
  ) +
  theme_minimal() +
  # Add this line to specify both limits and breaks for the x-axis
  scale_x_continuous(
    limits = c(0, 10),                 # Set the exact range of the x-axis
    breaks = seq(0, 10, by = 1)        # Set tick marks at every integer from 0 to 10
  ) +
  # Add this line to specify breaks
  theme(
    legend.position="none",
    legend.title = element_text(size = 8, face = "bold"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(size = 0.1),
    plot.title = element_text(hjust = 0.5, face='bold', size=12),
    axis.title=element_text(size=10,face="bold"),
    axis.text.y=element_text(size=8),
    axis.text.x=element_text(size=12),
    plot.margin = margin( 0.3,0.0,0.0,0.0, "cm")
  )




#Plot a grid the histograms of distances

list_distances_plt = list()
list_distances_plt[['Histogram_Actual_Data']] = plt_actual_distances
list_distances_plt[['Histogram_PermutedMedian_Data']] = plt_median_distances
n_cols_plot = 2
n_rows_plot = 1

width_page = 12
height_page = 6 * n_rows_plot #Inches
format = file_ext(output_file_name)

tmp_output_file_name = paste('Distances', output_file_name, sep = '.')

if(format == 'pdf' | format == 'PDF'){
  pdf(tmp_output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_distances_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page = plot_grid(plotlist = list_distances_plt, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
    ggsave(tmp_output_file_name, 
           plot = combined_page, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')






