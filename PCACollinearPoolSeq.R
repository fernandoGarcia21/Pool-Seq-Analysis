setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

#setwd('.')

#load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(FactoMineR)
library(cowplot)
library(foreach)
library(doParallel)
library(parallelly)
library("factoextra")

#library(FactoMineR, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")
#library(factoextra, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")

#library(factoextra, FactoMineR, lib.loc = "/nfs/scistore18/bartogrp/dgarciac/R/x86_64-pc-linux-gnu-library/4.3")

#install.packages("FactoMineR")

#install.packages("pacman")
#pacman checks whether the packages are missing and install them if needed
#pacman::p_load(ggplot2, FactoMineR, factoextra)


###########################################################
#                     Functions
###########################################################


######################################################
# Generate a scatter plot with icons for each group of snails
######################################################
plot_pca_figures <- function(res.pca, dataset, title){
  
  #Find the geographic location according to the sample name, to create the legend
  aux_n_cols = dim(dataset)[2]
  sample_names = data.frame(samples = colnames(dataset)[col_samples_start:aux_n_cols])
  sample_names$location <- sample_location_data$LOCATION[match(sample_names$samples, sample_location_data$SAMPLE)]
  #sample_names$Area <- sample_location_data$Area[match(sample_names$samples, sample_location_data$Sample)]
  variance_pc1 <- res.pca$eig[1,2]
  variance_pc2 <- res.pca$eig[2,2]
  
  
  #Create a dataframe with the coordinates of the principal components
  df_pca = as.data.frame(res.pca$ind$coord)
  df_pca$label = sample_names$location

  #Add the shape according to the population
  focus_pop_array = strsplit(focus_population, ',')[[1]]
  df_pca$shape_var <- ifelse(rownames(df_pca) %in% focus_pop_array, 17, 16) # 17 is a triangle, 16 is a circle. Change as needed.
  
  sahpes_locations = c('United Kingdom' = 9, 'Spain' = 8, 
                       'Ireland' = 3, 'Sweden' = 13, 
                       'Iceland' = 6, 'France' = 14, 
                       'Norway' = 1, 'Russia' = 4, 
                       'Faroe Islands' = 10, 
                       'Italy' = 18, 'Portugal' = 0,
                       'Austria' = 15, 'Korea' = 16)
  
  
  #Plot a scatterplot
  plt = ggplot(df_pca, aes(x = Dim.1, y = Dim.2*-1, 
                           shape = as.factor(shape_var),
                           color = as.factor(sample_names$location),
                           label = label)
                        )+
    geom_point(size = 3, alpha = 0.75)+
    #geom_text(hjust = 0.5, vjust = -0.5, size = 1.5, color = "red") + # Customize labels
    
    geom_text_repel(
      box.padding = 0.6, # Adjust padding around labels
      point.padding = 0.5, # Adjust padding around points
      segment.color = "black", # Color of connecting lines
      segment.size = 0.2, # Size of connecting lines
      min.segment.length = 0, # Ensure lines are always drawn
      arrow = arrow(length = unit(0.01, "npc")), #Add arrows if needed.
      force = 5, # Force the labels to move away from the points.
      size = 1.5,
      color = 'black',
      max.overlaps=length(df_pca$label)
    ) +
    
    theme_minimal()+
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size=10, face = "bold", hjust = 0.5),
      axis.title.x=element_text(size=8, face = "bold"),
      axis.title.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
    labs(x = paste('PC1 (', round(variance_pc1, digits=1), '%)', sep=""), 
         y = paste('PC2 (', round(variance_pc2, digits=1), '%)', sep=""))+
    ggtitle(title)+
    theme(legend.text = element_text(size = 5),
          legend.key.size = unit(0.2, "cm")) # Adjust size as needed
  
  return(plt)
}

######################################################
# Generate the PCA on a dataset
######################################################
generate_pca <- function(dataset){

  dataset_translated_transposed = t(dataset[,col_samples_start:dim(dataset)[2]])
  dataset_translated_transposed = matrix(as.numeric(dataset_translated_transposed), ncol = ncol(dataset_translated_transposed))
  snp_names = dataset[, 1]
  ind_names = colnames(dataset)[col_samples_start:dim(dataset)[2]]
  colnames(dataset_translated_transposed) <- snp_names
  rownames(dataset_translated_transposed) <- ind_names
  
  res.pca <- PCA(dataset_translated_transposed, graph = FALSE, ncp = 5)
  return(res.pca)
}


##################################################################
#Subsamples the snp dataset and generates a pca plot
##################################################################
create_pca_figure = function (p_replicate, p_genotypes_data){
  #To accelerate the testing, sample snp_sample_size SNPs randomly out of the total number of SNPs
  print(paste('Sampling', snp_sample_size, 'SNPs out of', dim(p_genotypes_data)[1]))
  tmp_random_snps = sample(dim(p_genotypes_data)[1], snp_sample_size, replace = FALSE)
  s_genotypes_data = p_genotypes_data[tmp_random_snps, ]
  

  #Invoke the function that performs a PCA on the given dataset
  collinear.res.pca <- generate_pca(s_genotypes_data)
  
  
  #Print the variance explained by the individuals in PC1 and PC2
  variance_contribution = as.data.frame(collinear.res.pca$ind$contrib[,c('Dim.1','Dim.2')],row.names = rownames(collinear.res.pca$ind$contrib))
  print('Sorted Variance explaned PC1:')
  print(variance_contribution[order(variance_contribution$Dim.1, decreasing = TRUE), ])
  
  print('Sorted Variance explaned PC2:')
  print(variance_contribution[order(variance_contribution$Dim.2, decreasing = TRUE), ])
  
  #Write two files with the PCA info for individuals
  #print('Write summary file with the PCA data')
  #aux_n_cols = dim(s_genotypes_data)[2]
  #sample_names = data.frame(samples = colnames(s_genotypes_data)[col_samples_start:aux_n_cols])
  #sample_names$location <- sample_location_data$Location[match(sample_names$samples, sample_location_data$Sample)]
  #sample_names$area <- sample_location_data$Area[match(sample_names$samples, sample_location_data$Sample)]
  #write.table(sample_names, file = paste0('Pools_', group_name, '_Replicate_', p_replicate, '.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  ind = get_pca_ind(collinear.res.pca)
  write.table(ind$coord, file = paste0('Pools_', group_name, '_PCs_Replicate_', p_replicate, '.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  write.table(ind$contrib, file = paste0('Pools_', group_name, '_Variance_Replicate_', p_replicate, '.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  tmp_title = paste('PoolSeq', group_name, 'Replicate No.', p_replicate)
  tmp_plt = plot_pca_figures(collinear.res.pca, s_genotypes_data, tmp_title)
  
  return(tmp_plt)
}


##################################################################
#Create a list of dot plots for the first n PCs after randomly smampling snps
##################################################################
create_multiple_pcs_plot = function (p_replicate, p_genotypes_data){
  
  list_pc_plots = list()
  #To accelerate the testing, sample 1000 SNPs randomly out of 10K
  print(paste('Sampling', snp_sample_size, 'SNPs out of', dim(p_genotypes_data)[1]))
  tmp_random_snps = sample(dim(p_genotypes_data)[1], snp_sample_size, replace = FALSE)
  s_genotypes_data = p_genotypes_data[tmp_random_snps, ]
  
  #s_genotypes_data = filter_maf(s_genotypes_data, col_samples_start)
  
  #Filter the snp dataset to remove incomplete samples (columns) and incomplete SNPs (rows)
  tmp_filter_by_areas = TRUE
  s_genotypes_data = remove_incomplete_samples_parallel(s_genotypes_data, max_na_snps, col_samples_start, log_file_name)
  s_genotypes_data = remove_incomplete_snps_parallel(s_genotypes_data, max_na_snails, max_na_snails_location, col_samples_start, sample_location_data, log_file_name)
  
  #Invoke the function that performs a PCA on the given dataset
  collinear.res.pca <- generate_pca(s_genotypes_data)
  
  res.pca = collinear.res.pca
  dataset = s_genotypes_data
  
  #Find the geographic location according to the sample name, to create the legend
  aux_n_cols = dim(dataset)[2]
  sample_names = data.frame(samples = colnames(dataset)[col_samples_start:aux_n_cols])
  sample_names$Location <- sample_location_data$Location[match(sample_names$samples, sample_location_data$Sample)]
  
  #Create a dataframe with three columns, the PC number, the PC value, and the location
  ind_pcs = as.data.frame(res.pca$ind$coord[,c('Dim.1','Dim.2','Dim.3','Dim.4')], 
                          row.names = rownames(res.pca$ind$contrib))
  ind_pcs$Location = sample_location_data$Location[match(rownames(ind_pcs), sample_location_data$Sample)]
  df_pcs = data.frame(PC = rep(c(1,2,3,4), each=nrow(ind_pcs)),
                      value = c(ind_pcs$Dim.1, ind_pcs$Dim.2, ind_pcs$Dim.3, ind_pcs$Dim.4), 
                      Location = rep(ind_pcs$Location, 4))
  
  #Create one dot plot for each of the first 4 PCs
  for(ip in 1:4){
    tmp_df_pc = df_pcs[df_pcs$PC == ip, ]
    
    #Remove points greather than or equal to the 95% quantile
    #tmp_df_pc = tmp_df_pc[tmp_df_pc$value < quantile(tmp_df_pc$value, 0.95) & tmp_df_pc$value > quantile(tmp_df_pc$value, 0.05), ]
    
    tmp_plt = ggplot(tmp_df_pc, aes(x=Location, y=value, colour = factor(Location))) +
      geom_jitter(shape=16, position=position_jitter(0.05), size = 2, alpha = 0.7) +
      theme_minimal()+
      ylab(paste('PC', ip, '- Replicate', p_replicate))+
      theme(legend.position = 'none',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    list_pc_plots[[paste0(p_replicate, '-', ip)]] = tmp_plt
  }
  
  return(list_pc_plots)
}


##################################################################
#Creates the PDF file with a grid of different PCAs, for the first two PCs
##################################################################
create_usual_pca_plot = function(p_genotypes_data){
  list_pca_plots = list() #Temporary store the PCA plots before plotting a grid
  
  #Iterate multiple replicates of the PCA
  for(i in 1:num_replicates){
    print(paste('Generating Replicate ',i))
    
    #Create a replicate PCA and add the plot to a list
    tmp_plot = create_pca_figure(i, p_genotypes_data)
    list_pca_plots[[i]] = tmp_plot
  }
  
  #Create grid with two plots per row
  plt = plot_grid(plotlist = list_pca_plots,
                  nrow = num_replicates,
                  ncol = 1
  )
  return (plt)
}


##################################################################
#Creates a grid with multiple replicates (rows) and the first n PCs
##################################################################
create_list_pcs_grid = function(p_genotypes_data){
  list_pca_plots = list() #Temporary store the PCA plots before plotting a grid
  
  #Iterate multiple replicates of the PCA
  for(i in 1:num_replicates){
    print(paste('Generating Replicate ',i))
    
    #Create a replicate PCA and add the plot to a list
    tmp_list_pc_plots = create_multiple_pcs_plot(i, p_genotypes_data)
    
    list_pca_plots = append(list_pca_plots, tmp_list_pc_plots)
  }
  
  #Create grid with two plots per row
  plt_grid = plot_grid(plotlist = list_pca_plots,
                       nrow = num_replicates,
                       ncol = 4
  )
  
  return(plt_grid)
  
}

###########################################################
#                Main body: Call to functions
###########################################################

#Input file paths
genotypes_file = 'Datasets/subsetEurope.frq' #SNP genotypes file
#genotypes_file = ''
sample_location_file = 'Datasets/Sample_Location.txt' #Two column file with the sample name and the geographical location
exclude_populations = c()
exclude_samples = c("17TFA",
                    "VA-3191-SPCAN-low_S22_L003",
                    "VA-3191-SPCAN-high_S23_L003",
                    "SPn_Crab_High",
                    "SPn_Wave_Low",
                    "SPs_Crab_High",
                    "SPs_Wave_Low",
                    "12POR")

num_replicates = 1 #Number of PCA replicates
group_name = 'Europe'
focus_population = 'VEN,3GAL'
snp_sample_size = 1000 #Number of SNPs to be sampled

#args = commandArgs(trailingOnly=TRUE)
## test if there is at least one argument: if not, return an error
#if (length(args) < 5) {
#  stop("At least five arguments must be supplied (input_file group focus_population focus_pool_population number_replicates sample_size_SNPs)", call.=FALSE)
#} else {
#  genotypes_file = args[1] #Frequencies estimated from sync file of popoolation
#  group_name = args[2] #Whether the group of the samples is Europe or North America
#  focus_population = args[3] #The name of the pool that will have a different shape (e.g. VEN)
#  num_replicates = as.numeric(args[4]) #How many PCA plots must be generated
#  snp_sample_size = as.numeric(args[5]) #How many SNPS must be sampled out of the total dataset
#}


#Read file and load data
genotypes_data_original = read.table(genotypes_file, header = TRUE, check.names = FALSE)
sample_location_data = read.table(sample_location_file, header = TRUE, sep = "\t", check.names = FALSE)

col_samples_start = 4 #The column number of the first sample

# Clear the log file before running (optional)
log_file_name = paste0('PCA.',group_name,'.log.txt')
writeLines("", log_file_name)


#Find the geographic location according to the sample name, to exclude unwanted populations and unwanted samples
aux_n_cols = dim(genotypes_data_original)[2]
sample_names = data.frame(samples = colnames(genotypes_data_original)[col_samples_start:aux_n_cols])
#sample_names$Location = sample_location_data$Location[match(sample_names$samples, sample_location_data$Sample)]

keep_samples = sample_names[!sample_names$samples %in% exclude_samples, 'samples'] #Find the rows (location names) that are not in the exclude list
genotypes_data = genotypes_data_original[,c(colnames(genotypes_data_original[1:2]),keep_samples)] #Filter the genotypes data to exclude unwanted locations


#Initialize figure params
output_file_name = ''
aux_width = 0
aux_height = 0


output_file_name = 'PoolSeq_PCA_2PCs_'
aux_width = 9
aux_height = 6*num_replicates
plt = create_usual_pca_plot(genotypes_data)


#Create a PDF file with the grid of PCA figures
#print('Creating PDF file')
output_file_name = paste0(output_file_name, '_', group_name, '_', snp_sample_size, '_', num_replicates, '.pdf')
#pdf(output_file_name, width = aux_width, height = aux_height)
plt
#dev.off()
print(paste('PCA completed:',output_file_name))


