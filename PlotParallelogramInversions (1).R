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
plot_scatter_parallelogram = function(p_df_frequencies, p_inv_name){
  p_title = p_inv_name
  
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


file_af_populations = 'Datasets/subsetNAmerica.frq';
file_af_dignostic_snps = 'Datasets/InvDiagnosticSNPs.minCount2.namerica_allLG.DP30.frq';
file_inv_coordinates = 'Datasets/Inversion coordinates Lsax_20250313.txt';

#Parallelogram approach to detect inversions on Pool-Seq data
populataions_set = c('SM-BOO-high-merged', 'SM-BOO-low-merged')


args = commandArgs(trailingOnly=TRUE)
# test if there is at least five arguments: if not, return an error
if (length(args) < 4) {
  stop("At least five arguments must be supplied (diagnostic_snps_file, frequencies_true_sax_population, frequencies_test_population, name_true_sax_population, name_test_sax_population)", call.=FALSE)
} else {
  file_af_populations = args[1] #Path of the file where the AF of the pools is found
  file_af_dignostic_snps = args[2] #Path of the file where the af of the diagnostic SNPs
  file_inv_coordinates = args[3] #Path of the file where the start and end coordinates of each inversion is found
  populations = args[4] #Name of the populations that will be analysed in pairs. Separated by comma.
}

print('Finish reading parameters')
populataions_set = strsplit(populations, ',')[[1]]

print('Reading the allele frequencies')
#Read the true saxatilis information
data_af = read.table(file_af_populations, header = TRUE, check.names = FALSE)

#Read the diagnostic SNPs af
data_diagnostic_snps_af = read.table(file_af_dignostic_snps, header = TRUE, check.names = FALSE)


data_coordinates = read.table(file_inv_coordinates, header = TRUE, check.names = FALSE)


#---------------------------------------------------------------------------------------------------------------------
#split the LG name in the af file and diagnostic SNPs file
data_af$LG = sapply(strsplit(data_af$chrom, "_"), function(x) x[1])
data_diagnostic_snps_af$LG = sapply(strsplit(data_diagnostic_snps_af$chrom, "_"), function(x) x[1])


#Split the last part of the inversion name by "." and get only the first part, which is the LGC
split_LG_name = sapply(strsplit(data_coordinates$inv, "\\."), function(x) x[1])
first_part = substr(split_LG_name, 1, 2)
remaining_part = substr(split_LG_name, 4, nchar(split_LG_name))
split_LG_name = paste0(first_part, remaining_part)
data_coordinates$LG = split_LG_name
data_coordinates$start = as.integer(data_coordinates$start)
data_coordinates$end = as.integer(data_coordinates$end)

print(colnames(data_coordinates))

#---------------------------------------------------------------------------------------------------------------------


#For each pair of population (rows) plot the AF within each coordinate window
pop_pairs = generate_unique_pairs(populataions_set);

#named list that will contain other lists, one list for each chromosome
listPlotsPopulations = list();

#For each pair of populations
for(tmp_pop_pair in pop_pairs){
  
  #Get the pair of populations
  tmp_cols_populations = c(tmp_pop_pair[1], tmp_pop_pair[2]);
  print(tmp_cols_populations)
  
  
  #named list that will contain all plots for the inversions of a population pair
  listPlotsInversions = list();
  cont_inversions = 0
  #Create one figure for each inversion
  for(i_inv in 1:dim(data_coordinates)[1]){
    
    tmp_inv = data_coordinates[i_inv,]
    tmp_af_inversion = data_af[data_af$LG == tmp_inv$LG 
                               & data_af$pos > tmp_inv$start-1
                               & data_af$pos < tmp_inv$end + 1, tmp_cols_populations]
    
    
    #if there are SNPs within the inversion
    if(nrow(tmp_af_inversion) > 4){
      #Since there are many SNPs within each invesion, use only a quarter of the SNP content
      tmp_sampled_rows = sample(dim(tmp_af_inversion)[1], dim(tmp_af_inversion)[1]/4, replace = FALSE)
      tmp_af_inversion = tmp_af_inversion[tmp_sampled_rows,]
      
      print(paste(tmp_inv$inv, 'with this number of SNPs:', dim(tmp_af_inversion)[1]))
      
      tmp_af_inversion$diagnostic = 'N'


      #Identify inversion diagnostic SNPs that belong to the inversion and add them to the dataset to plot in a different colour
      tmp_af_diagnostic_snps_inversion = data_diagnostic_snps_af[data_diagnostic_snps_af$LG == tmp_inv$LG & 
                                                                   data_diagnostic_snps_af$pos > tmp_inv$start-1 & 
                                                                   data_diagnostic_snps_af$pos < tmp_inv$end+1, tmp_cols_populations]
      
      if(dim(tmp_af_diagnostic_snps_inversion)[1] > 0){

        print(paste('Diagnostic SNPs with data ', dim(tmp_af_diagnostic_snps_inversion)[1]))
        
        tmp_af_diagnostic_snps_inversion$diagnostic = 'Y'
        tmp_af_inversion = rbind(tmp_af_inversion, tmp_af_diagnostic_snps_inversion)
        
      }#Close if diagnostic SNPs were found within the inversion coordinates
      
    }#Close if there are SNPs within the inversion coordinates
    else{
      # Add a third column with zero rows
      tmp_af_inversion$diagnostic <- character(0)
    }
    
    #Rename the last two columns of the dataset for convenience at the plotting function
    print(colnames(tmp_af_inversion))
    colnames(tmp_af_inversion) = c('popA', 'popB', 'Diagnostic')
    
    
    #---------------------------------------------------------------------------------------------------------------------
    
    
    #Generate a scatter plot of the frequencies
    tmp_plot = plot_scatter_parallelogram(tmp_af_inversion, tmp_inv$inv)
    
    cont_inversions = cont_inversions+1
    listPlotsInversions[[cont_inversions]] = tmp_plot
    
    
  } #Close for inversion coordinates
  
  #Add the list of chromosomes plots to the list of population pairs
  listPlotsPopulations[[paste(tmp_pop_pair[1], '-', tmp_pop_pair[2])]] = listPlotsInversions

}#close for population pairs


#Start plotting the inversions
resolution <- 300 # Adjust as needed (DPI)

for (page_name in names(listPlotsPopulations)) {
  page = listPlotsPopulations[[page_name]]
  
  n_cols_plot = 4 # Desired number of columns in the grid
  n_plots = length(page)
  
  if (n_plots > 0) {
    n_rows_plot = ceiling(n_plots / n_cols_plot) # Calculate rows correctly
    
    width_page = 2.5 * n_cols_plot #Inches
    height_page = 2.5 * n_rows_plot #Inches
    
    image_width = width_page * resolution #pixels
    image_height = height_page * resolution #pixels
    
    combined_page = plot_grid(plotlist = page, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")
    
    ggsave(paste0("Inversions_", page_name, ".jpg"), 
           plot = combined_page, 
           device = "jpeg", 
           units = "px",
           width = image_width, 
           height = image_height, 
           dpi = resolution)
  }
}


