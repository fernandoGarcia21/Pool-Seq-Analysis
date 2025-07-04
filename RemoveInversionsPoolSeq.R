setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')
#setwd('.')


args = commandArgs(trailingOnly=TRUE)

coordinates_file = 'Datasets/Inversion coordinates Lsax_20250313.txt' #Inversion scoring new txt file
input_file = ''


#Frequencies file, the output of the Popoolation pipeline
input_file = 'Datasets/InvDiagnosticSNPs.minCount2.namerica_allLG.DP30.frq'

# Get file name with extension
file_name_with_ext <- basename(input_file)

# Get file name without extension
file_name_without_ext <- tools::file_path_sans_ext(file_name_with_ext)

# Get file extension
file_ext <- tools::file_ext(file_name_with_ext)


input_data = read.table(input_file, header = TRUE)
coordinates = read.table(coordinates_file, header = TRUE)

#Transform the window value from text to numeric for practicity

coordinates$start = as.integer(coordinates$start)
coordinates$end = as.integer(coordinates$end)

#Split the last part of the inversion name by "." and get only the first part, which is the LGC
split_LG_name = sapply(strsplit(coordinates$inv, "\\."), function(x) x[1])
first_part = substr(split_LG_name, 1, 2)
remaining_part = substr(split_LG_name, 4, nchar(split_LG_name))
split_LG_name = paste0(first_part, remaining_part)
coordinates$LG = split_LG_name

#coordinates$nWindow = as.numeric(gsub(',','',coordinates$Window))

#Split the chromosome info for practicity
SNP_info = input_data[,c(1,2)]
SNP_info$LG = sapply(strsplit(SNP_info$chrom,"_"), `[`, 1)
SNP_info$index = 1:nrow(SNP_info)


#Go through every SNP and identify weather to keep it or to exclude it
#Exclude SNPs that fall within one inversion or it is uncertain
rows_keep_collinear = c()
rows_keep_inversions = c()

print(paste('Procesing', nrow(SNP_info) , 'SNPs'))



for(i_inv in 1:nrow(coordinates)){
  tmp_inv = coordinates[i_inv,]
  
  #Identify the indexes of the SNPs that follow within the current inversion in the loop
  indexes_snps_inversion_index = SNP_info[which(SNP_info$LG == tmp_inv$LG &
                                                  SNP_info$pos > tmp_inv$start -1 &
                                                  SNP_info$pos < tmp_inv$end +1), ]$index
  
  rows_keep_inversions = append(rows_keep_inversions, indexes_snps_inversion_index)
}

#The collinear Loci will be SNPs whose indexes were not found as inversion indexes
rows_keep_collinear = SNP_info$index
rows_keep_collinear = rows_keep_collinear[!rows_keep_collinear %in% rows_keep_inversions]



for(tmp_row_snp in 1:nrow(SNP_info)){
  #We just need the first two columns: LG and coordinate
  tmp_snp = SNP_info[tmp_row_snp, ]

  #Check if the snp falls within an inversion for the LG
  tmp_snp_inversion = coordinates[coordinates$LG == tmp_snp$LG &
                                    coordinates$start < tmp_snp$pos+1 &
                                    coordinates$end > tmp_snp$pos - 1, ]
  
  #If the window is not an inversion (No or '') the SNP is a in a collinear region
  if(nrow(tmp_snp_inversion) == 0){
    rows_keep_collinear = append(rows_keep_collinear, tmp_row_snp)
  }else{
      rows_keep_inversions = append(rows_keep_inversions, tmp_row_snp)
  }
}

print(paste('Initial SNPs:', dim(SNP_info)[1]))
print(paste('Discarded SNPs (Uncertain):', dim(SNP_info)[1] - (length(rows_keep_collinear) + length(rows_keep_inversions))))
print(paste('Collinear SNPs to write:', length(rows_keep_collinear)))
print(paste('Inversion SNPs to write:', length(rows_keep_inversions)))


#Write the output of Collinear loci to a file adding a filter name to the end of the file name
write.table(input_data[rows_keep_collinear, ], 
            file = paste0(file_name_without_ext,".collinear.", file_ext), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = '\t')

print(paste('The output file', paste0(file_name_without_ext,".collinear.", file_ext), 'was printed.'))

#Write the output of Inversions to a file adding a filter name to the end of the file name
write.table(input_data[rows_keep_inversions, ], 
            file = paste0(file_name_without_ext,".inversions.", file_ext), 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = '\t')

print(paste('The output file', paste0(file_name_without_ext,".inversions.", file_ext), 'was printed.'))

