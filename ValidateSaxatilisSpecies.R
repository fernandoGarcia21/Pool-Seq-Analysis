setwd('.')
setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(dplyr)
library(ggplot2)



######################################################
# Computes the FST for each SNP of two populations 
# and returns an array with FST values
######################################################
compute_fst <- function(aaf_pop_A, aaf_pop_B){
  raf_pop_A = 1 - aaf_pop_A
  raf_pop_B = 1 - aaf_pop_B
  
  expec_heterozygosity_A = 2 * raf_pop_A * aaf_pop_A
  expec_heterozygosity_B = 2 * raf_pop_B * aaf_pop_B
  
  h_s = (expec_heterozygosity_A + expec_heterozygosity_B) / 2
  h_t_ref = (raf_pop_A + raf_pop_B) / 2
  h_t_alt = (aaf_pop_A + aaf_pop_B) / 2
  h_t = 2 * h_t_ref * h_t_alt
  #Find the indexes of all ht equal 0
  h_t_0 = which(h_t == 0)
  
  fst = (h_t - h_s) / h_t
  print(h_t_0)
  #Replace those Nan values of fst with 0
  fst[h_t_0] = 0
  
  fst
}



#Usar los SNPs diagnosticos
file_species_diagnostic_snps = 'Datasets/SpeciesDiagnosticSNPs.txt';
file_af_populations = 'Datasets/subsetEurope.frq';
file_af_for_sampling = 'Datasets/subsetEurope.frq';

#Usar una población de la que estamos seguros es Saxatilis
true_sax_population = 'VA-3191-TSWANG-low_S24_L003' #Angklavebukten, Sweden
test_sax_population = 'VEN'#The SF populations RED and OAK

n_random_sampling = 500


args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 5) {
  stop("At least five arguments must be supplied (diagnostic_snps_file, frequencies_true_sax_population, frequencies_test_population, name_true_sax_population, name_test_sax_population)", call.=FALSE)
} else {
  file_species_diagnostic_snps = args[1] #List of diagnostic SNPs: two columns LG and position
  file_af_populations = args[2] #Path of the file where the AF of the truely saxatilis pool is found
  file_af_for_sampling = args[3] #Path of the file that I will use to sample random loci
  true_sax_population = args[4] #Name of the sample of the truely saxatilis population
  test_sax_population = args[5] #Name of the sample of the saxatilis population that we want to find out the species
}

print('Reading diagnostic SNPs')
data_diagnostic_snps = read.table(file_species_diagnostic_snps, header = FALSE, check.names = FALSE)
colnames(data_diagnostic_snps) = c("chrom", "pos")

print('Reading the allele frequencies')
#Read the true saxatilis information
data_af = read.table(file_af_populations, header = TRUE, check.names = FALSE)


#Extract only the data of the true saxatilis population and the diatnostic snps
data_true_sax_population = data_af[, c(colnames(data_af)[1:3], true_sax_population)]
filtered_true_sax_population = data_true_sax_population %>% inner_join(data_diagnostic_snps, by = c("chrom", "pos"))
print(paste('The SNPs in the filtered true dataset: ', dim(filtered_true_sax_population)[1]))


#Extract only the data of the test saxatilis population and the diatnostic snps
data_test_sax_population = data_af[, c(colnames(data_af)[1:3], test_sax_population)]
filtered_test_sax_population = data_test_sax_population %>% inner_join(data_diagnostic_snps, by = c("chrom", "pos"))
print(paste('The SNPs in the filtered test dataset: ', dim(filtered_test_sax_population)[1]))

#Eliminate the large dataset to release memmory
rm(data_af, data_true_sax_population, data_test_sax_population)

#Read the dataset that will be use to randomly sample SNPs
data_af_for_sampling = read.table(file_af_for_sampling, header = TRUE, check.names = FALSE)
# Sample without replacement
sampled_rows = sample(row_number(data_af_for_sampling), n_random_sampling, replace = FALSE)

#Obtain the frequencies of the sampled rows
tmp_sample_frequencies = data_af_for_sampling[sampled_rows, c(colnames(data_af_for_sampling)[1:3], c(true_sax_population, test_sax_population))]

#Remove the dataframe for sampling that we do not need anymore to release memmory
rm(data_af_for_sampling)



#Add two additional columns to estimate the sax allele, ie. the allele that is the most common in the Saxatilis population
filtered_true_sax_population$sax_allele = rep('', nrow(filtered_true_sax_population))
filtered_true_sax_population$freq_sax_allele = rep(0,nrow(filtered_true_sax_population))

filtered_test_sax_population$sax_allele = rep('', nrow(filtered_test_sax_population))
filtered_test_sax_population$freq_sax_allele = rep(0,nrow(filtered_test_sax_population))

#Identify the Sax allele (the most frequent in the true population)
for(i in 1:nrow(filtered_true_sax_population)){
  tmp_SNP = filtered_true_sax_population[i,]
  freq_sax_allele = tmp_SNP[4]
  name_sax_allele = 'R'
  
  #Get the default frequency of the reference allele for the test population
  freq_test_allele = filtered_test_sax_population[i,4]
  
  #If the current frequency is not the most common allele then the sax allele is the alternative
  if(freq_sax_allele < 0.5){
    freq_sax_allele = 1-freq_sax_allele
    name_sax_allele = 'A'
    
    freq_test_allele = 1-freq_test_allele
  }
  
  #Update the saxatilis allele in both populations
  filtered_true_sax_population[i, 'sax_allele'] = name_sax_allele;
  filtered_true_sax_population[i, 'freq_sax_allele'] = freq_sax_allele;
  
  filtered_test_sax_population[i, 'sax_allele'] = name_sax_allele;
  filtered_test_sax_population[i, 'freq_sax_allele'] = freq_test_allele;
}

#Join both frequencies in one DF
df_frequencies_populations = data.frame(`true_sax_population` = filtered_true_sax_population$freq_sax_allele,
                                        `test_sax_population` = filtered_test_sax_population$freq_sax_allele)

#Estimate the FST of the alternative allele
fst_estimates_alt_allele = compute_fst(filtered_true_sax_population[,4], filtered_test_sax_population[,4])


#Estimate the FST of the sampled loci
fst_estimates_sample_alt_allele = compute_fst(tmp_sample_frequencies[,4], tmp_sample_frequencies[,5])
quantiles_random_SNPs = quantile(fst_estimates_sample_alt_allele, c(0.025, 0.975))

df_fst_estimates = data.frame(SNP=paste0(filtered_true_sax_population$chrom,'_', filtered_true_sax_population$pos),
                              FST = fst_estimates_alt_allele,
                              mean_FST = mean(fst_estimates_sample_alt_allele),
                              upper_limit = as.numeric(quantiles_random_SNPs[1]),
                              lower_limit = as.numeric(quantiles_random_SNPs[2]))



################################################################
#         Plot FST and 95% quantile
################################################################

plt_FST = ggplot(df_fst_estimates, aes(x = SNP)) +
  geom_point(aes(y = FST)) + # Plot FST points
  geom_hline(aes(yintercept = mean_FST), color = "red", linewidth = 1) + # Mean FST line
  geom_hline(aes(yintercept = upper_limit), linetype = "dashed", color = "blue") + # Upper limit line
  geom_hline(aes(yintercept = lower_limit), linetype = "dashed", color = "blue") + # Lower limit line
  labs(
    x = "SNP",
    y = "FST",
    title = paste0('FST ', test_sax_population, ' vs ', true_sax_population)
  ) +
  theme_minimal()+ # Optional: Use a minimal theme for cleaner look
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # Rotate labels vertically

#Create a PDF file with the grid of PCA figures
print('Creating PDF file')
output_file_name = paste0('FST_', test_sax_population, '_', true_sax_population, '.pdf')
pdf(output_file_name, width = 6, height = 7)
plt_FST
dev.off()
print(paste('FST completed:',output_file_name))





################################################################
#         Plot Scatter Plot with correlation lines
################################################################
#Linear model of correlation
m_out = lm(true_sax_population~test_sax_population, df_frequencies_populations) #slope
r2_out = format(summary(m_out)$adj.r.squared, digits = 2) #R2 or the coefficient of determination

#Estimate the correlation coefficient and p value
corr <- cor.test(x=df_frequencies_populations$true_sax_population, y=df_frequencies_populations$test_sax_population, method = 'spearman')
corr

p_title = paste('Frequency Correlation of the Saxatilis Allele:', test_sax_population)

#Generate a ggplot object with thecorrelation between the Saxatilis allele frequencies of the populations
plt_correlation = ggplot(df_frequencies_populations, 
                      aes(x=true_sax_population, y=test_sax_population))+
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, linetype='F1', size=0.8)+
  annotate(geom="text",x=0.9,y=0.6,
           label=(paste0("R2==",r2_out)),
           parse=TRUE, size=3, fontface = 'bold')+
  ylim(-0.1, 1)+
  xlim(0, 1)+
  xlab(paste('Freq. Sax allele in ', true_sax_population))+
  ylab(paste('Freq. Sax allele in ', test_sax_population))+
  theme_minimal()+
  ggtitle(p_title)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10),
        #plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1))


#Create a PDF file with the grid of PCA figures
print('Creating PDF file')
output_file_name = paste0('Correlation_', test_sax_population, '.pdf')
pdf(output_file_name, width = 6, height = 5)
plt_correlation
dev.off()
print(paste('PCA completed:',output_file_name))

