setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)
library(tools)
library(RColorBrewer)
library(tidyr)

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

read_and_merge_files = function(p_file_list, p_sep){
  
  # Initialize an empty list to store the data frames
  data_frames_list <- list()
  # Loop through each file in the list
  for (file_path in p_file_list) {
    # Read the CSV file into a data frame
    current_df = read.table(file_path, header = TRUE, sep = p_sep, check.names = FALSE)
    
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
    cat("\nNo .fst-list.csv files found in the specified folder.\n")
    merged_df = NULL # Or however you want to handle the case with no files
  }
  return(merged_df)
}



#Directory where the diversity statistics are found
#fst_files_dir_path = 'Grenedalf/WindowFST Europe';
files_dir_fst_path = 'Grenedalf/WindowFST';
files_dir_diversity_path = 'Grenedalf/WindowDiversity';
files_dir_tajimasD_path = 'Grenedalf/TajimasD_100K';

col_start_data = 4

#output_file_name = 'Europe_PiRatioFST_Grenedalf.pdf';
output_file_name = 'Grenedalf_Distributions.jpg';

focus_population = 'VEN'
area_focus_population = list()
area_focus_population[['OAK']] = 'NA'
area_focus_population[['RED']] = 'NA'
area_focus_population[['VEN']] = 'EU'
populations_highlight = c('VEN','OAK','RED', 'EAS-low', 'SM-BOO-low-merged')

output_file_name = paste(focus_population, output_file_name, sep = '.')

# List all files in the specified folder with the .diversity.csv extension
file_list_fst = list.files(path = files_dir_fst_path, pattern = paste0("\\.",focus_population,".*\\.fst\\.csv$"), full.names = TRUE)
file_list_diversity = list.files(path = files_dir_diversity_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)
file_list_tajimasD = list.files(path = files_dir_tajimasD_path, pattern = paste0("\\.",area_focus_population[[focus_population]],".*\\.diversity\\.csv$"), full.names = TRUE)

merged_fst_df = read_and_merge_files(file_list_fst, '\t')
merged_diversity_df = read_and_merge_files(file_list_diversity, '\t')
merged_tajimasD_df = read_and_merge_files(file_list_tajimasD, '\t')




################################################################################
#        ******************* VIOLIN PLOTS *******************
################################################################################


################################################################################
#FST
################################################################################

tmp_plot_fst_df = merged_fst_df[,col_start_data:ncol(merged_fst_df)]
#Extract the names of the columns of the fst dataframe and update the name
tmp_colnames_fst = colnames(tmp_plot_fst_df)
tmp_colnames_fst = gsub("\\.fst$", "", tmp_colnames_fst) # Delete ".fst" at the end
tmp_colnames_fst = sapply(strsplit(tmp_colnames_fst, ":"), function(x) x[2])
colnames(tmp_plot_fst_df) = tmp_colnames_fst

#Flaten the dataframe
tmp_plot_fst_df <- pivot_longer(
  data = tmp_plot_fst_df,
  cols = everything(),
  names_to = "Population",
  values_to = "FST"
)

#Remove NA values
tmp_plot_fst_df = na.omit(tmp_plot_fst_df)

#Generate a list of intercalated colours to distinguish the chromosomes
n = length(tmp_colnames_fst)
color_array = generate_intercalated_colors(n, color1 = "#80cdc1", color2 = "#dfc27d")
sorted_pop_names = tmp_colnames_fst[order(tmp_colnames_fst)]
listColourCode = setNames(as.list(color_array), sorted_pop_names) 

# Create a new 'color' column by looking up the color for each population
tmp_plot_fst_df$Color = listColourCode[as.character(tmp_plot_fst_df$Population)]


#Estimate the medians
medians_fst_df = tmp_plot_fst_df %>% 
                  group_by(Population) %>% 
                  summarise(median = median(FST))


#Plot violins
plt_violin_fst = ggplot(tmp_plot_fst_df, aes(x = Population, y = FST, fill = Color, col=Color)) +
  geom_violin(alpha = 0.7) +
  geom_point(data = medians_fst_df, aes(x = Population, y = median), fill = '#262626', col = '#262626', alpha = 0.5)+
  labs(title = focus_population,
       x = "Population",
       y = "FST",
       fill = "Population") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face = "bold"),
        axis.title=element_text(size=8))+
  guides(fill = FALSE, col = FALSE)

#Print the plot
ggsave(paste0('Violin.FST.', output_file_name), 
       plot = plt_violin_fst, 
       device = "jpeg", 
       units = "in",
       width = 8, 
       height = 5, 
       dpi = 300)




################################################################################
#Pi: Diversity
################################################################################
list_colnames_pi = colnames(merged_diversity_df)
list_colnames_pi_number = grep("\\.theta_pi$", list_colnames_pi)
tmp_plot_pi_df = merged_diversity_df[,list_colnames_pi_number]

#Extract the names of the columns of the fst dataframe and update the name
tmp_colnames_pi = colnames(tmp_plot_pi_df)
tmp_colnames_pi = gsub("\\.theta_pi$", "", tmp_colnames_pi) # Delete ".theta_pi" at the end
colnames(tmp_plot_pi_df) = tmp_colnames_pi

#Flaten the dataframe
tmp_plot_pi_df <- pivot_longer(
  data = tmp_plot_pi_df,
  cols = everything(),
  names_to = "Population",
  values_to = "Pi"
)

#Remove NA values
tmp_plot_pi_df = na.omit(tmp_plot_pi_df)

#0.02 for Europe
tmp_plot_pi_df = tmp_plot_pi_df[tmp_plot_pi_df$Pi <= 0.05, ]
tmp_plot_pi_df = tmp_plot_pi_df[!is.na(tmp_plot_pi_df$Pi), ]

#Generate a list of intercalated colours to distinguish the chromosomes
n = length(tmp_colnames_pi)
color_array = generate_intercalated_colors(n, color1 = "#80cdc1", color2 = "#dfc27d")
sorted_pop_names = tmp_colnames_pi[order(tmp_colnames_pi)]
listColourCode = setNames(as.list(color_array), sorted_pop_names) 

# Create a new 'color' column by looking up the color for each population
tmp_plot_pi_df$Color = listColourCode[as.character(tmp_plot_pi_df$Population)]


#Estimate the medians
medians_pi_df = tmp_plot_pi_df %>% 
  group_by(Population) %>% 
  summarise(median = median(Pi))


#Plot violins
plt_violin_pi = ggplot(tmp_plot_pi_df, aes(x = Population, y = Pi, fill = Color, col=Color)) +
  geom_violin(alpha = 0.7) +
  geom_point(data = medians_pi_df, aes(x = Population, y = median), fill = '#262626', col = '#262626', alpha = 0.5)+
  labs(title = paste('Nucleotide Diversity', '-', area_focus_population[focus_population]),
       x = "Population",
       y = "Theta Pi",
       fill = "Population") +
  #ylim(0, max(tmp_plot_pi_df$Pi) + 0.05)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face = "bold"),
        axis.title=element_text(size=8))+
  guides(fill = FALSE, col = FALSE)

#Print the plot
ggsave(paste0('Violin.Pi.', output_file_name), 
       plot = plt_violin_pi, 
       device = "jpeg", 
       units = "in",
       width = 8, 
       height = 5, 
       dpi = 300)




################################################################################
#Tajima's D
################################################################################
list_colnames_tajimasD = colnames(merged_tajimasD_df)
list_colnames_tajimasD_number = grep("\\.tajimas_d$", list_colnames_tajimasD)
tmp_plot_tajimasD_df = merged_tajimasD_df[,list_colnames_tajimasD_number]

#Extract the names of the columns of the tajimas dataframe and update the name
tmp_colnames_tajimasD = colnames(tmp_plot_tajimasD_df)
tmp_colnames_tajimasD = gsub("\\.tajimas_d$", "", tmp_colnames_tajimasD) # Delete ".tajimas_d" at the end
colnames(tmp_plot_tajimasD_df) = tmp_colnames_tajimasD

#Flaten the dataframe
tmp_plot_tajimasD_df <- pivot_longer(
  data = tmp_plot_tajimasD_df,
  cols = everything(),
  names_to = "Population",
  values_to = "TajimasD"
)

#Remove NA values
tmp_plot_tajimasD_df = na.omit(tmp_plot_tajimasD_df)

#Generate a list of intercalated colours to distinguish the chromosomes
n = length(tmp_colnames_tajimasD)
color_array = generate_intercalated_colors(n, color1 = "#80cdc1", color2 = "#dfc27d")
sorted_pop_names = tmp_colnames_tajimasD[order(tmp_colnames_tajimasD)]
listColourCode = setNames(as.list(color_array), sorted_pop_names) 

# Create a new 'color' column by looking up the color for each population
tmp_plot_tajimasD_df$Color = listColourCode[as.character(tmp_plot_tajimasD_df$Population)]

#Estimate the medians
medians_tajimasD_df = tmp_plot_tajimasD_df %>% 
  group_by(Population) %>% 
  summarise(median = median(TajimasD))

#Plot violins
plt_violin_tajimasD = ggplot(tmp_plot_tajimasD_df, aes(x = Population, y = TajimasD, fill = Color, col=Color)) +
  geom_violin(alpha = 0.7) +
  geom_point(data = medians_tajimasD_df, aes(x = Population, y = median), fill = '#262626', col = '#262626', alpha = 0.5)+
  labs(title = paste("Tajima's D", '-', area_focus_population[focus_population]),
       x = "Population",
       y = "Tajima's D",
       fill = "Population") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=8, face = "bold"),
        axis.title=element_text(size=8))+
  guides(fill = FALSE, col = FALSE)

#Print the plot
ggsave(paste0('Violin.TajimasD.', output_file_name), 
       plot = plt_violin_tajimasD, 
       device = "jpeg", 
       units = "in",
       width = 8, 
       height = 5, 
       dpi = 300)









################################################################################
#        ******************* HISTOGRAMS *******************
################################################################################

#Lists of plots
list_fst_plots = list()
list_pi_plots = list()
list_tajimasD_plots = list()

################################################################################
#FST
################################################################################
list_populations_fst = colnames(merged_fst_df)
fst_pairs = gsub("\\.fst$", "", list_populations_fst) # Delete ".fst" at the end
fst_pairs = gsub(":", " vs ", fst_pairs)       # Replace ":" with "vs"

for(i in (col_start_data : length(fst_pairs))){
  tmp_title = fst_pairs[i]
  tmp_fst_data = merged_fst_df[,i]
  #colnames(tmp_fst_data) = c('fst')
  
  fill_color = "#85b1ff"
  if(tmp_title %in% paste(focus_population, 'vs', populations_highlight)){
    fill_color = "#f30057"
  }
  
  tmp_plt_fst = ggplot(data.frame(fst = tmp_fst_data), aes(x = fst)) +
    geom_histogram(fill = fill_color, color = "black") +
    xlim(-0.1,1.1)+
    ylim(0,1800)+
    labs(title = tmp_title,
         x = "FST",
         y = "Frequency") +
    theme_minimal()+
    theme(plot.margin = unit(c(0.5, 0, 0.2, 0), "cm"),# Top, Right, Bottom, Left
          plot.title = element_text(hjust = 0.5, size=8,),
          axis.title=element_text(size=6)) 
  
  list_fst_plots[[tmp_title]] = tmp_plt_fst
}


################################################################################
#PI
################################################################################
list_populations_pi = colnames(merged_diversity_df)
list_populations_pi_number = grep("\\.theta_pi$", list_populations_pi)


for(i in list_populations_pi_number){
  
  tmp_title = gsub("\\.theta_pi$", "", list_populations_pi[i])
  tmp_pi_data = merged_diversity_df[,i]
  
  fill_color = "#dfc27d"
  if(tmp_title %in% populations_highlight){
    fill_color = "#f30057"
  }
  
  tmp_plt_pi = ggplot(data.frame(pi = tmp_pi_data), aes(x = pi)) +
    geom_histogram(binwidth = 0.01, fill = fill_color, color = "black") +
    xlim(0,0.4)+
    ylim(0,2000)+
    labs(title = tmp_title,
         x = "Theta Pi",
         y = "Frequency") +
    theme_minimal()+
    theme(plot.margin = unit(c(0.5, 0, 0.2, 0), "cm"),# Top, Right, Bottom, Left
          plot.title = element_text(hjust = 0.5, size=8,),
          axis.title=element_text(size=6)) 
  
  list_pi_plots[[tmp_title]] = tmp_plt_pi
}




################################################################################
#Tajimas's D
################################################################################
list_populations_tajimasD = colnames(merged_tajimasD_df)
list_populations_tajimasD_number = grep("\\.tajimas_d$", list_populations_tajimasD)


for(i in list_populations_tajimasD_number){
  
  tmp_title = gsub("\\.tajimas_d$", "", list_populations_tajimasD[i])
  tmp_tajimasD_data = merged_tajimasD_df[,i]
  
  fill_color = "#80cdc1"
  if(tmp_title %in% populations_highlight){
    fill_color = "#f30057"
  }
  
  tmp_plt_tajimasD = ggplot(data.frame(tajimasD = tmp_tajimasD_data), aes(x = tajimasD)) +
    geom_histogram(binwidth = 0.1, fill = fill_color, color = "black") +
    xlim(-4,1)+
    ylim(0,12500)+
    labs(title = tmp_title,
         x = "Tajima's D",
         y = "Frequency") +
    theme_minimal()+
    theme(plot.margin = unit(c(0.5, 0, 0.2, 0), "cm"),# Top, Right, Bottom, Left
          plot.title = element_text(hjust = 0.5, size=8,),
          axis.title=element_text(size=6)) 
  
  list_tajimasD_plots[[tmp_title]] = tmp_plt_tajimasD
}




#Generate a PDF with the distributions

n_cols_plot = 8
n_rows_plot = ceiling(length(list_fst_plots)/n_cols_plot)

width_page = n_cols_plot * 3
height_page = 3 * n_rows_plot #Inches
format = file_ext(output_file_name)

if(format == 'pdf' | format == 'PDF'){
  pdf(output_file_name, width = width_page, height = height_page)
  plot_grid(plotlist = list_fst_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
  dev.off()
}else{
  if(format == 'jpg' | format == 'JPG' | format == 'JPEG' | format == 'jpeg'){
    
    combined_page_fst = plot_grid(plotlist = list_fst_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    ggsave(paste0('FST.', output_file_name), 
           plot = combined_page_fst, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
    
    combined_page_pi = plot_grid(plotlist = list_pi_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    ggsave(paste0('Pi.', output_file_name), 
           plot = combined_page_pi, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
    
    combined_page_tajimasD = plot_grid(plotlist = list_tajimasD_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    ggsave(paste0('TajimasD.', output_file_name), 
           plot = combined_page_tajimasD, 
           device = "jpeg", 
           units = "in",
           width = width_page, 
           height = height_page, 
           dpi = 300)
  }
  
  
}

print('Output file created')



