setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(dplyr)



####################################################################################################

# Additional filter to remove outliers from the inversions

####################################################################################################
identify_outlier_indices_sd_grouped <- function(data_group, sd_threshold = 2) {
  mean_val <- mean(data_group$Frequency)
  sd_val <- sd(data_group$Frequency)
  upper_bound <- mean_val + sd_threshold * sd_val
  lower_bound <- mean_val - sd_threshold * sd_val
  outlier_indices_within_group <- which(data_group$Frequency > upper_bound | data_group$Frequency < lower_bound)
  # Return the original row indices within the dataframe for these outliers
  original_indices <- data_group$Index[outlier_indices_within_group]
  return(as.numeric(original_indices))
}


filter_af_by_outliers = function (p_data_diagnostic_snps_af, p_populations) {
  
  #Identify outliers of frequencies by population, inversion and arrangement
  df_outlier_freqs = data.frame()
  for(tmp_population in p_populations){
    tmp_af_population = p_data_diagnostic_snps_af[,c(tmp_population, 'Inversion', 'pos', 'Arrangement')]
    colnames(tmp_af_population) = c('Frequency','Inversion', 'Position', 'Arrangement')
    
    tmp_af_population$Index = rownames(tmp_af_population)
    
    # Apply the function grouped by Inversion and Arrangement and get the indexes
    # of frequencies that are beyond 4 standard deviations of the data within the arrangement
    outlier_indices_grouped_sd <- tmp_af_population %>%
      group_by(Inversion, Arrangement) %>%
      summarise(Outlier_Indices = list(identify_outlier_indices_sd_grouped(cur_data(), 4))) %>%
      ungroup()
    
    
    tmp_outliers_df = tmp_af_population[unlist(outlier_indices_grouped_sd$Outlier_Indices),]
    tmp_outliers_df$Population = tmp_population
    df_outlier_freqs = rbind(df_outlier_freqs, tmp_outliers_df)
  }
  
  #Exclude the outliers from the df datast
  unique_index_outliers = unique(df_outlier_freqs$Index)
  # Create a logical vector of row numbers to keep
  rows_to_keep_logical <- !(rownames(p_data_diagnostic_snps_af) %in% unique_index_outliers)
  
  final_data_diagnostic_snps_af = p_data_diagnostic_snps_af[rows_to_keep_logical,]
  return(final_data_diagnostic_snps_af)
}



#Estimate and write inversion frequencies for all populations
file_af_dignostic_snps = 'Datasets/InvDiagnosticSNPs.minCount2.europe_allLG.200.DP30_MAF002.frq';
file_output_clean_afs = 'Datasets/InvDiagnosticSNPs.minCount2.europe_allLG.200.DP30_MAF002.clean.frq';

file_inv_coordinates = 'Datasets/Inversion coordinates Lsax_20250313.txt';
file_snp_arrangements = 'Datasets/parallelogram.new.trimmed.200.coordinates.arr.txt';


#Read the diagnostic SNPs af
data_diagnostic_snps_af = read.table(file_af_dignostic_snps, header = TRUE, check.names = FALSE)

#Read the coordinates of the inversion breakpoints
data_coordinates = read.table(file_inv_coordinates, header = TRUE, check.names = FALSE)

#Read the SNP-Arrangement associations
data_snp_arrangements = read.table(file_snp_arrangements, header = FALSE, check.names = FALSE)
colnames(data_snp_arrangements) = c('chrom','pos','Arrangement')

col_start_freq = 4
populations = colnames(data_diagnostic_snps_af)[col_start_freq:dim(data_diagnostic_snps_af)[2]]


#Split the last part of the inversion name by "." and get only the first part, which is the LGC
split_LG_name = sapply(strsplit(data_coordinates$inv, "\\."), function(x) x[1])
first_part = substr(split_LG_name, 1, 2)
remaining_part = substr(split_LG_name, 4, nchar(split_LG_name))
split_LG_name = paste0(first_part, remaining_part)
data_coordinates$LG = split_LG_name
data_coordinates$start = as.integer(data_coordinates$start)
data_coordinates$end = as.integer(data_coordinates$end)

#Add the LG to the diagnostic SNPs
data_diagnostic_snps_af$LG = sapply(strsplit(data_diagnostic_snps_af$chrom, "_"), function(x) x[1])

#Add the arrangement information to the frequencies dataset
data_diagnostic_snps_af = left_join(data_diagnostic_snps_af, data_snp_arrangements, by = c("chrom", "pos"))

#Add the inversion information to the frequencies dataset
data_diagnostic_snps_af$Inversion = NA


#Search for the inversion that corresponds to each diagnostic SNP position
for(i_diag in 1:dim(data_diagnostic_snps_af)[1]){
  tmp_snp = data_diagnostic_snps_af[i_diag,]
  
  tmp_inv = data_coordinates[data_coordinates$LG == tmp_snp$LG
                             & data_coordinates$start <= tmp_snp$pos
                             & data_coordinates$end >= tmp_snp$pos, ]

  if(dim(tmp_inv)[1] < 1){
    print('SNP outside inversions')
    print(tmp_snp[,1:2])
  }else{
    data_diagnostic_snps_af[i_diag, ]$Inversion = tmp_inv[1,]$inv
  }
}



#Remove diagnostic SNPs that were not allocated in an inversion
data_diagnostic_snps_af = data_diagnostic_snps_af[!is.na(data_diagnostic_snps_af$Inversion),]

#Exclude the outliers in the distributions of arrangement frequencies
clean_data_diagnostic_snps_af = filter_af_by_outliers(data_diagnostic_snps_af, populations)

##########################################################################################
#Re-write the inversion diagnostic frequencies
write.table(file=file_output_clean_afs, clean_data_diagnostic_snps_af[,1:(ncol(clean_data_diagnostic_snps_af)-3)], quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
##########################################################################################






#Estimate the mean of the inversions for each population
df_mean_freqs = data.frame(Inversion = c(), Arrangement = c(), Frequency=c(), Population=c())
for(tmp_population in populations){
  tmp_af_population = clean_data_diagnostic_snps_af[,c(tmp_population, 'Inversion', 'Arrangement')]
  colnames(tmp_af_population) = c('Population','Inversion','Arrangement')
  
  means_diagnostics <- as.data.frame(tmp_af_population %>%
                                       group_by(Inversion, Arrangement) %>%
                                       summarize(Frequency = mean(Population)))
  means_diagnostics$Population = tmp_population
  df_mean_freqs = rbind(df_mean_freqs, means_diagnostics)
}



#Write the mean frequencies in a file
write.table(file='Datasets/InversionFrequencies.Namerica.clean.txt', df_mean_freqs, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


