setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

# Install Bioconductor if you haven't already
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("GenomicRanges", "IRanges", "regioneR"))

library(GenomicRanges)
library(IRanges)
library(regioneR) # For permutation tests
library(ggplot2)
library(dplyr)


#Focus population
focus_population = 'VEN'
array_min_locations = c(2,3,4,5)

# --- 1. Load Outlier Windows ---
outlier_windows_df <- read.table(paste0(focus_population, ".LowPiNoThreshold.sharedFstLowPiRegionsOutlierWindows.txt"), header = TRUE, sep = "\t")
outlier_windows_df$chr = sapply(strsplit(outlier_windows_df$chrom, "_"), function(x) x[1])


# Convert to GRanges object
outlier_windows_gr <- GRanges(
  seqnames = outlier_windows_df$chr,
  ranges = IRanges(start = outlier_windows_df$start, end = outlier_windows_df$end)
)


# --- 2. Load all SNPs from the high-low shore project ---
all_snps_df <- read.table("Datasets/AnjaOutliers/reduced.POOLHL14_OUTLIERS_invs_filt_full.txt", header = TRUE, sep = " ")
colnames(all_snps_df) = c('cp','outlier_count')

all_snps_df$chr = sapply(strsplit(all_snps_df$cp, "_"), function(x) x[1])
all_snps_df$pos = as.numeric(sapply(strsplit(all_snps_df$cp, "_"), function(x) x[2]))

#keep only collinear loci by excluding NA's from pos
all_snps_df = na.omit(all_snps_df)


all_snps_gr <- GRanges(
  seqnames = all_snps_df$chr,
  ranges = IRanges(start = all_snps_df$pos, end = all_snps_df$pos)
)




################################## ONLY LG WITH CANDIDATE WINDOWS *****************************

total_significant_snps_df = data.frame()
total_non_significant_snps_df = data.frame()

#Identify significant enrichment for different numbers of locations where a SNP is high-low outlier
for(min_loc_number in array_min_locations){
  
  #Subset only the outlier SNPs (outliers in n or more locations)
  outlier_snps_df <- all_snps_df[all_snps_df$outlier_count > min_loc_number-1, ]
  
  
  # Convert to GRanges object (SNPs are point regions, so start = end)
  outlier_snps_gr <- GRanges(
    seqnames = outlier_snps_df$chr,
    ranges = IRanges(start = outlier_snps_df$pos, end = outlier_snps_df$pos),
    snp_id = outlier_snps_df$cp
  )
  
  #Temporal df to store the p values
  results_pvalue_df = data.frame()
  
  seqlevels(outlier_windows_gr)
  tmp_lg_list = intersect(seqlevels(outlier_windows_gr), seqlevels(outlier_snps_gr))
  
  
  for(LG in tmp_lg_list){
    #Exclude SNPs that are not part of the LGs of candidate regions under selection
    all_snps_LG_gr <- keepSeqlevels(all_snps_gr, LG, pruning.mode="coarse")
    outlier_snps_LG_gr <- keepSeqlevels(outlier_snps_gr, LG, pruning.mode="coarse")
    
    
    # Calculate overlaps
    overlaps_snps_in_windows <- findOverlaps(outlier_snps_LG_gr, outlier_windows_gr)
    a <- length(unique(queryHits(overlaps_snps_in_windows))) # Outlier SNPs in outlier windows
    
    # Outlier SNPs NOT in outlier windows
    b <- length(outlier_snps_LG_gr) - a
    
    # Other SNPs in outlier windows
    other_snps_gr <- subsetByOverlaps(all_snps_LG_gr, outlier_snps_LG_gr, invert = TRUE)
    overlaps_other_snps_in_windows <- findOverlaps(other_snps_gr, outlier_windows_gr)
    c <- length(unique(queryHits(overlaps_other_snps_in_windows)))
    
    # Other SNPs NOT in outlier windows
    d <- length(other_snps_gr) - c
    
    
    # Create contingency table
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                                dimnames = list(c("Outlier SNP", "Other SNP"),
                                                c("In Outlier Window", "Not In Outlier Window")))
    
    
    #print(contingency_table)
    
    # Perform Fisher's Exact Test
    fisher_test_result <- fisher.test(contingency_table, alternative = "greater")
    cat(LG, ': ', fisher_test_result$p.value, '\n')
    # A small p-value indicates significant enrichment (odds ratio > 1).
    
    results_pvalue_df = rbind(results_pvalue_df, c(LG, a, b, c, d, fisher_test_result$p.value))
    
    #If the p-value is significant, store the outlier SNPs to be printed in the end
    if(fisher_test_result$p.value <= 0.05){
      # Get the GRanges object of outlier SNPs that overlap with outlier windows
      if(length(overlaps_snps_in_windows) > 0){
        outlier_snps_in_outlier_windows_gr <- outlier_snps_LG_gr[unique(queryHits(overlaps_snps_in_windows))]
        significant_snps_df <- as.data.frame(outlier_snps_in_outlier_windows_gr)
        significant_snps_df$num_locations = min_loc_number
        
        total_significant_snps_df = rbind(total_significant_snps_df, significant_snps_df)
      }
    }else{
      # Get the GRanges object of outlier SNPs that overlap with outlier windows but the overlapping is not statistically signifcant
      if(length(overlaps_snps_in_windows) > 0){
        outlier_snps_in_outlier_windows_gr <- outlier_snps_LG_gr[unique(queryHits(overlaps_snps_in_windows))]
        non_significant_snps_df <- as.data.frame(outlier_snps_in_outlier_windows_gr)
        non_significant_snps_df$num_locations = min_loc_number
        
        total_non_significant_snps_df = rbind(total_non_significant_snps_df, non_significant_snps_df)
      }
    }
    
    
    
    
  }#Close for LGs
  
  colnames(results_pvalue_df) = c('LG',"a In Outlier Window", "b Not In Outlier Window","c In Outlier Window", "d Not In Outlier Window",'p-value')
  results_pvalue_df$`p-value` = as.numeric(results_pvalue_df$`p-value`)
  results_pvalue_df = results_pvalue_df[order(results_pvalue_df$`p-value`),]
  
  results_pvalue_df$`p-value (formatted)` <- formatC(
    results_pvalue_df$`p-value`,
    format = "e",   # 'e' for scientific notation
    digits = 3      # Number of decimal places
  )
  
  #Write the p-values of enrichment for a given number of locations where a snp is outlier
  write.table(results_pvalue_df, 
              file = paste0(focus_population,".minLoc.", min_loc_number, ".Fisher.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  
}# Close for location numbers


#Finally write the list of significantly enriched SNPs

# Filter for unique snp_id keeping the maximum num_locations
#filtered_significant_snps_df <- total_significant_snps_df[,c('snp_id','num_locations')] %>%
#  group_by(snp_id) %>%       # Group the dataframe by snp_id
#  filter(num_locations == min(num_locations)) %>% # For each group, keep only rows where num_locations is the min
#  ungroup()                  # Ungroup the dataframe (good practice after grouping operations)




write.table(total_significant_snps_df[,c('snp_id','num_locations')], 
            file = paste0(focus_population, ".significantSNPs.Fisher.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(total_non_significant_snps_df[,c('snp_id','num_locations')], 
            file = paste0(focus_population, ".nonSignificantSNPs.Fisher.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


####################### CALCULATE FISHER TEST FOR ALL LG TOGER


## Calculate overlaps
#overlaps_snps_in_windows <- findOverlaps(outlier_snps_gr, outlier_windows_gr)
#a <- length(unique(queryHits(overlaps_snps_in_windows))) # Outlier SNPs in outlier windows
#
## Outlier SNPs NOT in outlier windows
#b <- length(outlier_snps_gr) - a
#
## Other SNPs in outlier windows
#other_snps_gr <- subsetByOverlaps(all_snps_gr, outlier_snps_gr, invert = TRUE)
#overlaps_other_snps_in_windows <- findOverlaps(other_snps_gr, outlier_windows_gr)
#c <- length(unique(queryHits(overlaps_other_snps_in_windows)))
#
## Other SNPs NOT in outlier windows
#d <- length(other_snps_gr) - c
#
#
## Create contingency table
#contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
#                            dimnames = list(c("Outlier SNP", "Other SNP"),
#                                            c("In Outlier Window", "Not In Outlier Window")))
#print(contingency_table)
#
## Perform Fisher's Exact Test
#fisher_test_result <- fisher.test(contingency_table, alternative = "greater")
#print(fisher_test_result)
## A small p-value indicates significant enrichment (odds ratio > 1).
