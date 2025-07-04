#setwd('.')
setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(dplyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(Cairo)

#Generate non-repeted pairs of elements from an array
generate_unique_pairs = function(items) {
  n = length(items)
  
  if (n < 2) {
    return(list()) # Return an empty list if there are fewer than 2 items
  }
  
  pairs = list()
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pairs = append(pairs, list(c(items[i], items[j])))
    }
  }
  
  return(pairs)
}

#Generates a Scatter plot with the allele frequencies of two populations
plot_scatter_parallelogram = function(p_df_frequencies, p_window_start, p_window_end){
  p_title = paste(p_window_start, '-', p_window_end)
  
  #Generate a ggplot with the allele frequencies of the populations
  plt_parallelogram = ggplot(p_df_frequencies, 
                           aes(x=popA, y=popB, color = Diagnostic))+
    geom_point(size = 0.2, alpha = 0.3, shape = 1)+
    ylim(-0.1, 1)+
    xlim(0, 1)+
    theme_minimal()+
    ggtitle(p_title)+
    theme(
          #plot.margin = unit(c(0,0.5,0.5,0), "cm"),
          plot.title = element_text(size = rel(0.5)),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.1),
          legend.position = "none")+
    scale_colour_manual(values = c('N'='gray', 'Y'='red'))+
    guides(color = guide_legend(show = FALSE))
  
  return(plt_parallelogram)
}

# Function to arrange plots in a row with row name
arrange_row_with_name <- function(row_name, plots, p_max_number_windows) {
  perc_label_width = 0.025
  perc_figures_width = ((1-perc_label_width) / p_max_number_windows) * length(plots)
  print(paste('Figures Width',perc_figures_width))
  
  if (length(plots) == 0) {
    return(NULL)
  }
  
  row_label <- ggdraw() + draw_label(row_name, x = 0, hjust = 0, vjust = 0.5, size = 12) #Create label
  arranged_plots <- plot_grid(plotlist = plots, ncol = length(plots)) #Arrange plots
  
  cols_to_complete = p_max_number_windows - length(plots)
  rel_width_empty_col = ((1-perc_label_width) / p_max_number_windows) * cols_to_complete
  if(cols_to_complete > 0){
    row_empty_cell <- ggdraw()
    plot_grid(row_label, arranged_plots, row_empty_cell,  ncol = 3, rel_widths = c(perc_label_width, perc_figures_width, rel_width_empty_col)) #Combine label and plots.
    
  }else{
    plot_grid(row_label, arranged_plots, ncol = 2, rel_widths = c(perc_label_width, perc_figures_width)) #Combine label and plots.
  }
  
}





file_af_populations = 'Datasets/subsetNAmerica.frq';
file_inv_diagnostic_snps = 'Datasets/TranslatedSNPs.txt';

#Parallelogram approach to detect inversions on Pool-Seq data
slidingWindow = 2000000 #1 million base pairs
populataions_set = c('SM-BOO-high-merged', 'SM-BOO-low-merged')


# args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
# if (length(args) < 3) {
#   stop("At least five arguments must be supplied (diagnostic_snps_file, frequencies_true_sax_population, frequencies_test_population, name_true_sax_population, name_test_sax_population)", call.=FALSE)
# } else {
#   file_af_populations = args[1] #Path of the file where the AF of the pools is found
#   slidingWindow = as.numeric(args[2]) #Size of the window of SNP positions
#   populations = args[3] #Name of the populations that will be analysed in pairs. Separated by comma.
#   populataions_set = strsplit(populations, ',')[[1]]
# }

print('Finish reading parameters')


step_size = slidingWindow+1

print('Reading the allele frequencies')
#Read the true saxatilis information
data_af = read.table(file_af_populations, header = TRUE, check.names = FALSE)
data_diagnostic_snps = read.table(file_inv_diagnostic_snps, header = TRUE, check.names = FALSE)
data_diagnostic_snps = data_diagnostic_snps[!is.na(data_diagnostic_snps$Translated_pos), ]

#---------------------------------------------------------------------------------------------------------------------
#split the LG name in the af file
data_af$LG = sapply(strsplit(data_af$chrom, "_"), function(x) x[1])
#---------------------------------------------------------------------------------------------------------------------


#Get the min and max positions for each chromosome in the frequencies dataset
data_chr <- data_af %>%
  group_by(chrom) %>%
  summarise(min_position = min(pos),
            max_position = max(pos))

#Exclude chromosomes where there are not enough snps
data_chr = data_chr[data_chr$max_position > data_chr$min_position,]

#Generate the chunks of positions for a given slidingWindow within each chromosome
resultWindows <- data_chr %>%
  rowwise() %>%
  mutate(
    windows = list(seq(min_position, max_position - slidingWindow, by = step_size))
  ) %>%
  tidyr::unnest(windows) %>%
  mutate(window_end = windows + slidingWindow) %>%
  select(chrom, windows, window_end)

resultWindows = as.data.frame(resultWindows)


#For each pair of population (rows) plot the AF within each coordinate window
pop_pairs = generate_unique_pairs(populataions_set);

#named list that will contain other lists, one list for each chromosome
listPlotsPopulations = list();
max_number_windows = resultWindows %>% 
                      group_by(chrom) %>%
                      summarise(nWindows = n()) %>% 
                      slice_max(nWindows, n = 1)

#For each pair of populations
for(tmp_pop_pair in pop_pairs){
  #For each chromosome
  tmp_cols_populations = c(tmp_pop_pair[1], tmp_pop_pair[2]);
  print(tmp_cols_populations)
  
  #named list that will contain all plots for the windows within each chromosome, every element of the list is a chromosome
  listPlotsChromosomes = list();
  
 for(tmp_chr in data_chr$chrom){
   print(tmp_chr)
   tmp_chr_windows = resultWindows[resultWindows$chrom == tmp_chr, ];
   #For each window
   
   #named list that will contain all plots for the windows of a chromosome
   listPlotsWindows = list();
   
   cont_windows = 0
   for(tmp_i_window in seq(1:nrow(tmp_chr_windows))){
    tmp_window = tmp_chr_windows[tmp_i_window, ]

    print(tmp_window)
    
    #Subset the frequencies dataset to get only the SNPs in the window for just the pair of pulations
    tmp_data_pops_window = data_af[data_af$chrom == tmp_window$chrom & 
                                     data_af$pos > tmp_window$windows & 
                                     data_af$pos < tmp_window$window_end,
                                   tmp_cols_populations
                                   ]
    
    
    
    #Since there are many SNPs within each simdow, use only a quarter of the SNP content
    tmp_sampled_rows = sample(dim(tmp_data_pops_window)[1], dim(tmp_data_pops_window)[1]/4, replace = FALSE)
    tmp_data_pops_window = tmp_data_pops_window[tmp_sampled_rows,]
    
    
    #---------------------------------------------------------------------------------------------------------------------
    
    #Only if there are SNPs within the downsampled window, add the diagnostic SNPs
    if(dim(tmp_data_pops_window)[1] > 0){
      tmp_data_pops_window$diagnostic = 'N'
      #print(dim(tmp_data_pops_window))
      
      
      #split the LG name in the window 
      tmp_window$LG = sapply(strsplit(tmp_window$chrom, "_"), function(x) x[1])
      
      #Identify inversion diagnostic SNPs that belong to the window and add them to the dataset to plot in a different colour
      tmp_diagnostic_snps_window = data_diagnostic_snps[data_diagnostic_snps$LG == tmp_window$LG & 
                                                          data_diagnostic_snps$Translated_pos > tmp_window$windows & 
                                                          data_diagnostic_snps$Translated_pos < tmp_window$window_end, ]
      
      if(dim(tmp_diagnostic_snps_window)[1] > 0){
        print(paste('Diagnostic SNPs found in Window ', dim(tmp_diagnostic_snps_window)[1]))
        
        tmp_data_inv_diagnostic_snps = data_af[data_af$LG %in% tmp_diagnostic_snps_window$LG 
                                                & data_af$pos %in% tmp_diagnostic_snps_window$Translated_pos, 
                                                tmp_cols_populations]
        
        print(paste('Diagnostic SNPs with data ', dim(tmp_data_inv_diagnostic_snps)[1]))
        
        #If diagnostic SNPs have AF data, add them to the plot dataframe
        if(dim(tmp_data_inv_diagnostic_snps)[1] > 0){
          tmp_data_inv_diagnostic_snps$diagnostic = 'Y'
          tmp_data_pops_window = rbind(tmp_data_pops_window, tmp_data_inv_diagnostic_snps)
        }
        
      }
      
     
    }else{
      # Add a third column with zero rows
      tmp_data_pops_window$diagnostic <- character(0)
    }
    
    #Rename the last two columns of the dataset for convenience at the plotting function
    colnames(tmp_data_pops_window) = c('popA', 'popB', 'Diagnostic')
   

    #---------------------------------------------------------------------------------------------------------------------
    
    
    #Generate a scatter plot of the frequencies
    tmp_plot = plot_scatter_parallelogram(tmp_data_pops_window, tmp_window$windows, tmp_window$window_end)
    
    cont_windows = cont_windows+1
    listPlotsWindows[[cont_windows]] = tmp_plot
    
   }#Close for windows
   
   #Add the list of windows plots of a chromosome to the list of chromosomes
   listPlotsChromosomes[[tmp_chr]] = listPlotsWindows
 }#Close for chromosomes
  
  #Add the list of chromosomes plots to the list of population pairs
  listPlotsPopulations[[paste(tmp_pop_pair[1], '-', tmp_pop_pair[2])]] = listPlotsChromosomes
}#Close for populations



resolution <- 300      # Adjust as needed (DPI)
# Generate the PDF
for (page_name in names(listPlotsPopulations)[1]) {
  
  #Obtain the data within each pair of populations
  page = listPlotsPopulations[[page_name]]
  
  #Estimate the height and width of the page according to the number of chromosomes within the page and the number of maximum windows
  max_num_plots = as.numeric(max_number_windows[2])
  height_page = 2*length(names(page))
  width_page = 2.5*max_num_plots
  
  max_dim = (50*resolution) - 100
  
  
  image_width <- max_dim  # Adjust as needed (pixels)
  image_height <- 400*length(names(page)) # Adjust as needed (pixels)
  
  #Generate a Pdf image
  #pdf(paste0("Scatter_", page_name, ".pdf"), width = width_page, height = height_page) # Adjust width and height as needed
  
  # Generate the PNG image
  #png(filename = paste0("Scatter_", page_name, ".png"),width = image_width, height = image_height, res = resolution)
  
  #CairoPNG(paste0("Scatter_", page_name, ".png"), width = image_width, height = image_height, res = resolution)
  
  #Get the rows with the figures per chromosome
  rows <- lapply(names(page), function(row_name) arrange_row_with_name(row_name, page[[row_name]], max_num_plots)) #Apply function to each row.
  
  # Combine rows into a single page
  if (length(rows) > 0) {
    combined_page <- plot_grid(plotlist = rows, ncol = 1, align = "v")
    ggsave(paste0("Scatter_", page_name, ".jpg"), plot = combined_page, device = "jpeg", units = "px",
           width = image_width, height = image_height, dpi = resolution)
    #print(combined_page)
  }
  
  #dev.off() # Close the PDF device
  
}



