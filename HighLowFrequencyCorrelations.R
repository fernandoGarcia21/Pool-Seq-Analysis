setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(ggplot2)
library(dplyr)
library(cowplot)


#Focus population
focus_population = 'VEN'
array_min_locations = c(2,3,4,5)

introduced_freq_list = list()
introduced_freq_list[['OAK']] = 'SanFrancisco'
introduced_freq_list[['RED']] = 'SanFrancisco'
introduced_freq_list[['VEN']] = 'VEN'

# --- 1.1. Load Introduced population frequencies ---
introduced_frequencies_df <- read.table(paste0("Datasets/minCount2.",introduced_freq_list[[focus_population]],".Collinear.frq"), header = TRUE, sep = " ")
colnames(introduced_frequencies_df) = gsub(".FREQ", "", colnames(introduced_frequencies_df))
introduced_frequencies_df$lg = sapply(strsplit(introduced_frequencies_df$chrom, "_"), function(x) x[1])
introduced_frequencies_df$cp = paste(introduced_frequencies_df$lg, introduced_frequencies_df$pos, sep= '_')

# --- 1.2. Load H-L Outliers frequencies ---
hl_outlier_frequencies_df <- read.table("Datasets/AnjaOutliers/olmc.2.POOLHL14_OUTLIERS_invs_filt_full.txt", header = TRUE, sep = "")

# --- 2. Load H-L outliers enriched in candidate selected regions  ---
hl_outlier_in_windows_df <- read.table(paste0(focus_population, ".nonSignificantSNPs.Fisher.txt"), header = TRUE, sep = "\t")

# --- 3. Create a dataframe with the relation between population, shore level, frequencies, and outlier status  ---
colnames_hl_freq = colnames(hl_outlier_frequencies_df)
col_names_hl_out =  grep("_out$", colnames_hl_freq, value = TRUE)
df_freq_columns = data.frame(out = col_names_hl_out, freq_high = paste0(substr(col_names_hl_out, 1, 4),'H'), freq_low = paste0(substr(col_names_hl_out, 1, 4),'L'))

# --- 4. Iterate over the min count outliers to estimate the frequency of the high shore allele  ---
list_plots = list()
for (tmp_min_loc_count in array_min_locations){
  tmp_hl_outlier_in_windows_mc_df = hl_outlier_in_windows_df[hl_outlier_in_windows_df$num_locations == tmp_min_loc_count,]
  if(nrow(tmp_hl_outlier_in_windows_mc_df) > 0){
    # subset the frequencies only for SNPs in candidate regions under selection
    tmp_hl_frequencies_mc_df = hl_outlier_frequencies_df[hl_outlier_frequencies_df$cp %in% tmp_hl_outlier_in_windows_mc_df$snp_id, ]
    
    
    #Filter the introduced population frequencies
    tmp_introduced_outliers = introduced_frequencies_df[introduced_frequencies_df$cp %in% tmp_hl_outlier_in_windows_mc_df$snp_id, ]
    
    
    #Create a dataframe with the SNP position, the high allele info, and the frequency of the high allele
    high_shore_alleles_df = data.frame()
    
    if(nrow(tmp_introduced_outliers) > 0){
    
      #Iterate over all SNPs in the high-shore frequencies dataset for a given set of minimum location counts
      for (i in seq(1:nrow(tmp_hl_frequencies_mc_df))){
        tmp_high_allele = 'R'
        tmp_avg_high_allele_freq = 0
        tmp_avg_low_allele_freq = 0
        
        tmp_snp_frequencies = tmp_hl_frequencies_mc_df[i,]
        
        #Idenfity the columns of the high shore populations that are outliers
        tmp_cols_freq_high = df_freq_columns$freq_high[which(t(tmp_snp_frequencies[,df_freq_columns$out]) ==TRUE)]
        tmp_frequencies_high = unlist(tmp_snp_frequencies[, tmp_cols_freq_high]) #extract the frequencies of the reference allele in the outlier populations high shore
        
        #Idenfity the columns of the low shore populations that are outliers
        tmp_cols_freq_low = df_freq_columns$freq_low[which(t(tmp_snp_frequencies[,df_freq_columns$out]) ==TRUE)]
        tmp_frequencies_low = unlist(tmp_snp_frequencies[, tmp_cols_freq_low]) #extract the frequencies of the reference allele in the outlier populations low shore
        
        
        #Use only outliers where the high allele is the same in all populations
        if(all(tmp_frequencies_high >= tmp_frequencies_low) | all(!(tmp_frequencies_high >= tmp_frequencies_low))){
          #Estimate the mean frequencies of reference allele in high and low shores
          tmp_avg_high_allele_freq = mean(tmp_frequencies_high) # estimate the average frequency of the ref allele
          tmp_avg_low_allele_freq = mean(tmp_frequencies_low) # estimate the average frequency of the ref allele
          
          #get the SNP from the introduced population frequencies
          tmp_introduced_snp = tmp_introduced_outliers[tmp_introduced_outliers$cp == tmp_snp_frequencies$cp, ]
          #print(tmp_introduced_snp)
          tmp_introduced_frequency = NA
          if(nrow(tmp_introduced_snp) > 0){
            tmp_introduced_frequency = tmp_introduced_snp[1,focus_population]
          }
          
          #Identify the allele that is more common in high than in low shore
          if(tmp_avg_low_allele_freq >= tmp_avg_high_allele_freq){
            tmp_high_allele = 'A'
            tmp_avg_high_allele_freq = 1 - tmp_avg_high_allele_freq
            
            if(!is.na(tmp_introduced_frequency)){
              tmp_introduced_frequency = 1 - tmp_introduced_frequency
            }
          }
          
          #Temporarly store the high-shore frequencies and snp information
          tmp_snp_info = strsplit(tmp_snp_frequencies$cp, "_")
          high_shore_alleles_df = rbind(high_shore_alleles_df, data.frame(cp = tmp_snp_frequencies$cp,
                                                                          lg = sapply(tmp_snp_info, function(x) x[1]),
                                                                          pos = sapply(tmp_snp_info, function(x) x[2]),
                                                                          h_allele = tmp_high_allele,
                                                                          h_freq = tmp_avg_high_allele_freq,
                                                                          i_h_freq = tmp_introduced_frequency))
        }else{
          print(paste('SNP', tmp_snp_frequencies$cp, 'excluded because the high-shore allele is not the same in all populations'))
        }#close else check same allele in all populations
      }#Close for all snps enriched in candidate regions under selection
      
    
    
      #Exclude NA values, for example SNPs that were not found in the frequencies of the introduced populations
      high_shore_alleles_df = na.exclude(high_shore_alleles_df)
      
      #write favored values
      write.table(high_shore_alleles_df, 
                  file = paste0(focus_population,".favoredAlleles.", tmp_min_loc_count, ".txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      high_shore_alleles_df_sorted = high_shore_alleles_df[order(high_shore_alleles_df$h_freq),]
      high_shore_alleles_df_sorted$index = seq(1:nrow(high_shore_alleles_df_sorted))
      
      
      plot_data_frame = high_shore_alleles_df_sorted[,c('index','h_freq')]
      colnames(plot_data_frame) = c('Locus','Frequency')
      
      tmp_df_introduced = high_shore_alleles_df_sorted[,c('index','i_h_freq')]
      colnames(tmp_df_introduced) = c('Locus','Frequency')
      
      plot_data_frame = rbind(plot_data_frame, tmp_df_introduced)
      plot_data_frame$type = rep(c('High-shore', focus_population), each  = nrow(high_shore_alleles_df_sorted))
      
      
      # Create a named vector for your colors
      # Start with known colors
      custom_colors <- c("High-shore" = "#80cdc1")
      
      # Add the color for focus_population using its value as the name
      custom_colors[focus_population] <- "#f30057"
      
      #Create plot
      tmp_plt = ggplot(plot_data_frame, aes(x = Locus, y = Frequency, fill = type)) +
        geom_point(
          alpha = 0.7,      # Transparency of the fill color
          shape = 21,       # Use a shape that has a fill and a border (e.g., 21-25)
          size = 3,         # Adjust point size if needed
          stroke = 0.2,       # Thickness of the border line
          color = "black"   # Color of the border (this 'color' applies to the border,
          # while the 'color' in aes() applies to the fill of shape 21)
        ) + # Add scatter points, alpha for transparency
        labs(
          title = paste0("High-shore allele (", tmp_min_loc_count, ' or more locations)'),
          x = "Locus",
          y = "Frequency",
          color = "Location" # Legend title
        ) +
        ylim(0,1)+
        theme_minimal() + # A clean theme
        scale_color_manual(values =custom_colors)+ # Custom colors
        theme(legend.position = 'bottom',
              plot.title = element_text(hjust = 0.5, face='bold', size = 12),
              legend.title = element_blank()
        )
      
      list_plots[[paste('Plt',tmp_min_loc_count)]] = tmp_plt
      
    }#Close if SNPs from high-low outliers were found in the pool-seq dataset of the introduced populations 
    
  }
}


n_cols_plot = length(array_min_locations)
n_rows_plot = 1

width_page = 4*n_cols_plot
height_page = 4 #Inches
output_file_name = paste0(focus_population,'.HighShoreFrequencyCorrelation.jpg')

combined_page = plot_grid(plotlist = list_plots, ncol = n_cols_plot, nrow = n_rows_plot, align = "v")

ggsave(output_file_name, 
       plot = combined_page, 
       device = "jpeg", 
       units = "in",
       width = width_page, 
       height = height_page, 
       dpi = 300)


print('Output file created')







