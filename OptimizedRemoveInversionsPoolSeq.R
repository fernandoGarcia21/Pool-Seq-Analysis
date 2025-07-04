setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(data.table)
library(tools)

coordinates_file <- 'Datasets/Inversion coordinates Lsax_20250313.txt'
input_file <- 'Datasets/InvDiagnosticSNPs.minCount2.namerica_allLG.DP30.frq'

# Get file name with extension
file_name_with_ext <- basename(input_file)
file_name_without_ext <- file_path_sans_ext(file_name_with_ext)
file_ext <- tools::file_ext(file_name_with_ext)

# Read data as data.tables
input_data <- fread(input_file)
coordinates <- fread(coordinates_file)

# Preprocess coordinates
coordinates[, start := as.integer(start)]
coordinates[, end := as.integer(end)]

split_LG_name <- sapply(strsplit(coordinates$inv, "\\."), function(x) x[1])
first_part <- substr(split_LG_name, 1, 2)
remaining_part <- substr(split_LG_name, 4, nchar(split_LG_name))
coordinates[, LG := paste0(first_part, remaining_part)]

# Preprocess SNP data
SNP_info <- input_data[, .(chrom, pos)]
SNP_info[, LG := sapply(strsplit(chrom, "_"), `[`, 1)]

# Convert to data.tables for efficient operations
setDT(SNP_info)
setDT(coordinates)

# Create a key for fast joins
setkey(coordinates, LG, start, end)

# Perform non-equi join to find overlaps
inversion_overlaps <- foverlaps(
  SNP_info,
  coordinates,
  by.x = c("LG", "pos", "pos"),
  by.y = c("LG", "start", "end"),
  type = "any",
  nomatch = 0L
)

# Identify SNPs in inversions
inversion_rows <- inversion_overlaps[, unique(i.chrom_pos)]
SNP_info[, chrom_pos := paste(chrom, pos, sep = "_")]
inversion_rows_index <- which(SNP_info$chrom_pos %in% inversion_rows)
collinear_rows_index <- which(!(SNP_info$chrom_pos %in% inversion_rows))

# Output
print(paste('Initial SNPs:', nrow(SNP_info)))
print(paste('Discarded SNPs (Uncertain):', nrow(SNP_info) - (length(collinear_rows_index) + length(inversion_rows_index))))
print(paste('Collinear SNPs to write:', length(collinear_rows_index)))
print(paste('Inversion SNPs to write:', length(inversion_rows_index)))

# Write output files
fwrite(input_data[collinear_rows_index], file = paste0(file_name_without_ext, ".collinear.", file_ext), sep = '\t', quote = FALSE)
print(paste('The output file', paste0(file_name_without_ext, ".collinear.", file_ext), 'was printed.'))

fwrite(input_data[inversion_rows_index], file = paste0(file_name_without_ext, ".inversions.", file_ext), sep = '\t', quote = FALSE)
print(paste('The output file', paste0(file_name_without_ext, ".inversions.", file_ext), 'was printed.'))