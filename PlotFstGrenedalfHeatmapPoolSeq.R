setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

#install.packages('pheatmap')
#install.packages("pheatmap", lib = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3/")
library("pheatmap")
library(ggplot2)
library(reshape2) # For melting the matrix
library(dplyr)
library(grid) # To explicitly control page breaks
library(RColorBrewer)


#Estimate and write inversion frequencies for all populations
fst_files_dir_path = 'Grenedalf/WholeFST Europe';
#fst_files_dir_path = 'Grenedalf/WholeFST Namerica';

#output_file_name = 'Europe_FST_GrenedalfHeatmap.pdf';
output_file_name = 'NorthAmerica_FST_GrenedalfHeatmap.pdf';


# List all files in the specified folder with the .fst-list.csv extension
file_list = list.files(path = fst_files_dir_path, pattern = "\\.fst-list\\.csv$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_frames_list <- list()

# Loop through each file in the list
for (file_path in file_list) {
  # Read the CSV file into a data frame
  current_df = read.table(file_path, header = TRUE, sep = '\t', check.names = FALSE)
  
  # Add the current data frame to the list
  data_frames_list[[length(data_frames_list) + 1]] = current_df
  
  # Optional: Print the name of the file that was just loaded
  cat("Loaded file:", basename(file_path), "\n")
}

# Merge all data frames in the list into a single data frame
if (length(data_frames_list) > 0) {
  merged_df = do.call(rbind, data_frames_list)
  cat("\nSuccessfully merged", length(data_frames_list), "files into a single data frame.\n")
  
  colnames(merged_df) = c("Population1", "Population2", "FST")
} else {
  cat("\nNo .fst-list.csv files found in the specified folder.\n")
  merged_df = NULL # Or however you want to handle the case with no files
}


#Estimate the average FST (groupping by poopulation1 and population1)
fst_df = merged_df %>% group_by(Population1, Population2) %>%
  summarise(FST = mean(FST))





#################################################################
#Plot the heatmap
#################################################################

# 1. Convert the FST dataframe into a FST matrix

#Append the complementary parts of the matrix to make it symetrical
tmp_complement = fst_df[,c(2,1,3)]
colnames(tmp_complement) = c("Population1", "Population2", "FST")
fst_melted = rbind(fst_df, tmp_complement)

fst_matrix = with(fst_melted, tapply(FST, list(Population1, Population2), FUN = identity))

# Determine the order of populations based on FST
# We need to find a way to order the populations based on their overall FST values
# in the lower triangle. One approach is to calculate the average FST
# for each population across the comparisons in the lower triangle.

# Get unique population names
all_populations <- unique(c(as.character(fst_melted$Population2), as.character(fst_melted$Population2)))

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
  coord_equal()+
  geom_text(aes(label = sprintf("%.2f", FST)), size = 3) # Add text layer



#Create a heatmap that includes a dendrogram using the euclidean distances of the FST estimats
plt_heatmapTree = pheatmap(fst_matrix, 
                           fontsize_row = 6,
                           fontsize_col = 6,
                           cutree_rows = 4, 
                           cutree_cols = 4, 
                           display_numbers = TRUE,
                           number_color = 'black',
                           color = brewer.pal(n = 11, name = "BrBG"))



ggsave('heatmapEurope.jpg', 
       plot = plt_heatmapTree, 
       device = "jpeg", 
       units = "in",
       width = 12, 
       height = 11, 
       dpi = 300)

pdf(output_file_name, width = 9, height = 8)
print(plt_heatmap)
# Force a new page before plotting pheatmap
grid.newpage()
print(plt_heatmapTree)



#Create a barplot with the FST of Venice or San Francisco vs the other populations
focus_populations = c('VEN','RED','OAK')
fst_focus_pops = fst_df[fst_df$Population1 %in% focus_populations, c(1,2,3)]
tmp_second_part = fst_df[fst_df$Population2 %in% focus_populations, c(2,1,3)]
names(tmp_second_part) = names(fst_focus_pops)
fst_focus_pops = rbind(fst_focus_pops, tmp_second_part)
unique_focus_pops = unique(fst_focus_pops$Population1)


for(tmp_pop in unique_focus_pops){
  
  fst_focus_pop = fst_focus_pops[fst_focus_pops$Population1 == tmp_pop,]
  fst_focus_pop = fst_focus_pop[order(fst_focus_pop$FST, decreasing = TRUE), ]
  
  fst_focus_pop$Population2 <- factor(fst_focus_pop$Population2, levels = unique(fst_focus_pop$Population2))
  
  plt_barplot = ggplot(fst_focus_pop, aes(x = Population2, y = FST, fill = FST)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = brewer.pal(n = 11, name = "BrBG"), name = "FST") + # Using a fixed number of colors for the gradient
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(), # Hide the x-axis label
          plot.title = element_text(hjust = 0.5) # Center the title (optional)
    )+
    labs(title = tmp_pop) # Add the plot title
  
  print(plt_barplot)
}

dev.off()

print('Output file created')





