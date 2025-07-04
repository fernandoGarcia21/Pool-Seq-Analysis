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



#Directory where the ourliers of the FST and viersity diversity statistics are found
files_dir_diversity_outliers_path = 'Grenedalf/WindowOutliers_FST_PI_10K';

#Directory where the FST windows are found, we just need this to know how many windows were used
files_dir_fst_path = 'Grenedalf/WindowFST_10K';
files_dir_diversity_path = 'Grenedalf/WindowDiversity_10K';

#output_file_name = 'Europe_PiRatioFST_Grenedalf.pdf';
output_file_name = 'SharedWindowsSanFrancisco.EastCoast.jpg';

#focus_population = 'Maine'
focus_population = 'East'
area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['East']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['PLY-high']] = 'NA'
area_focus_population[['Maine']] = 'NA'
area_focus_population[['SM-BOO-low-merged']] = 'NA'
area_focus_population[['VEN']] = 'EU'
area_focus_population[['EAS-low']] = 'EU'

output_file_name = paste(focus_population, output_file_name, sep = '.')


#Read outliers file OAK.SM-BOO-low-merged.WindowOutliers
file_list_diversity_outliers = list.files(path = files_dir_diversity_outliers_path, pattern = paste0("^",focus_population,".*\\.LowPiNoThreshold\\.sharedFstLowPiRegionsOutlierWindows\\.txt$"), full.names = TRUE)
if(length(file_list_diversity_outliers) > 0){
  diversity_outliers_df = read_and_merge_files(file_list_diversity_outliers, '\t')
}else{
  diversity_outliers_df = data.frame()
}


# List all files in the specified folder with the .diversity.csv extension
file_list_fst = list.files(path = files_dir_fst_path, pattern = paste0("\\.",focus_population,".*\\.fst\\.csv$"), full.names = TRUE)
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)

merged_fst_df = read_and_merge_files(file_list_fst, '\t')
merged_diversity_df = read_and_merge_files(file_list_diversity, '\t')

# Merge the dataframes by the first three columns
merged_df_all_windows = merge(merged_fst_df, merged_diversity_df, by = c("chrom", "start", "end"))
merged_df_all_windows = merged_df_all_windows[,1:3]
# Create a unique identifier for all windows DF
merged_df_all_windows <- merged_df_all_windows %>%
  mutate(unique_id = paste(chrom, start, end, sep = "_"))

#identify Unique comparison pairs
tmp_unique_groups = unique(diversity_outliers_df$group)
# Add new columns with FALSE values and add TRUE to the windows that are true in each group of outliers
for (group in tmp_unique_groups) {
  merged_df_all_windows[[group]] = FALSE
  tmp_outliers_group = diversity_outliers_df[diversity_outliers_df$group == group,]
  merged_df_all_windows[merged_df_all_windows$unique_id %in% tmp_outliers_group$unique_id, group] = TRUE
}


# Create a unique identifier
diversity_outliers_df <- diversity_outliers_df %>%
  mutate(unique_id = paste(chrom, start, end, sep = "_"))

# Function to create the list for ggvenn
create_venn_list <- function(group_name) {
  diversity_outliers_df %>%
    filter(group == group_name) %>%
    pull(unique_id)
}

# Get unique groups
groups <- unique(diversity_outliers_df$group)

# Create the list
venn_list <- lapply(groups, create_venn_list)
names(venn_list) <- groups # Name the list elements

# Print the list (for debugging)
print(venn_list)

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


########################################################################
# Print the shared outliers and estimate the probability of being shared
#########################################################################
if(length(shared_outliers) > 0){
  
  shared_diversity_outliers_df = diversity_outliers_df[diversity_outliers_df$unique_id %in% shared_outliers, c('chrom','start','end')]
  shared_diversity_outliers_df = unique(shared_diversity_outliers_df)
  #Export the shared outlier windows
  tmp_shared_outlier_file_name = paste(focus_population, 'sharedTopOutlierWindows.txt', sep = '.')
  write.table(shared_diversity_outliers_df, file = tmp_shared_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  
  #Estimate the probability of finding those shared outliers by chance
  #length(shared_outliers)
  #lengths_outliers = lengths(venn_list)
  #tmp_probabilities = length(shared_outliers) / lengths_outliers
  #probability_shared = prod(tmp_probabilities)
  #print(paste('The probability of sharing', length(shared_outliers), 'outliers is', probability_shared*100))
  
  
  num_outliers_group = diversity_outliers_df %>% 
    group_by(group) %>%
    summarise(Count = dplyr::n())
  
  
  observed_shared_outliers = length(shared_outliers) # The number of shared outliers you observed
  number_groups = nrow(num_outliers_group)
  group_cols = num_outliers_group$group
  count_chance_shared_occurrences = 0
  perm_replicates = 10
  
  for(i in 1:perm_replicates){
    tmp_shuffled_windows = merged_df_all_windows
    
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


#Get Unique windows for each population
library(purrr)
 unique_elements_purrr <- imap(venn_list, function(current_elements, current_group_name) {
   other_elements <- unlist(venn_list[names(venn_list) != current_group_name])
   setdiff(current_elements, other_elements)
 })
 print("\nElements unique to each group (using purrr):")
 print(unique_elements_purrr)
 
 for (tmp_unique_vs_pop in names(unique_elements_purrr)){
   tmp_unique_elements_pop = unique_elements_purrr[[tmp_unique_vs_pop]]
   
   unique_diversity_outliers_df = diversity_outliers_df[diversity_outliers_df$unique_id %in% tmp_unique_elements_pop, c('chrom','start','end')]
   #Export the shared outlier windows
   tmp_unique_outlier_file_name = paste(tmp_unique_vs_pop, focus_population, 'uniqueOutlierWindows.txt', sep = '.')
   write.table(unique_diversity_outliers_df, file = tmp_unique_outlier_file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
   
   
 }
 
