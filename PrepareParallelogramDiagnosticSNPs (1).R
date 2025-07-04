setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

library(tools)

# Function to process text files in a directory
process_text_files <- function(directory, output_file) {
  # Get a list of all .txt files in the directory
  files <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
  
  # Check if any .txt files were found
  if (length(files) == 0) {
    stop("No .txt files found in the specified directory.")
  }
  
  # Initialize an empty data frame to store the merged data
  merged_data <- data.frame()
  
  # Read and merge the content of each file
  for (file in files) {
    tryCatch({
      current_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if(dim(current_data)[1] > 0){
        
        file_name_no_extension <- file_path_sans_ext(file)
        
        # Split the file name by "_" and extract the last part which is the inversion name
        split_values_inv_name = sapply(strsplit(file_name_no_extension, "_"), function(x) x[4])
        
        #Split the last part of the file name by "." and get only the first part, which is the
        split_LG_name = sapply(strsplit(split_values_inv_name, "\\."), function(x) x[1])
        first_part = substr(split_LG_name, 1, 2)
        remaining_part = substr(split_LG_name, 4, nchar(split_LG_name))
        split_LG_name = paste0(first_part, remaining_part)
        
        current_data$Inv_Name = split_values_inv_name
        current_data$LG = split_LG_name
        colnames(current_data) = c('SNP','Inv_Name','split_LG_name')
        
        merged_data <- rbind(merged_data, current_data)
      
      }#Only if the file has data otherwise print the file is empty
      else{
        print(paste('The file', file, 'is empty'))
      }
    }, error = function(e) {
      warning(paste("Error reading file:", file, "-", e$message))
    })
  }
  
  # Check if merged_data is empty after reading all files.
  if (nrow(merged_data) == 0){
    stop("No data was read from the text files.")
  }
  
  # Split the first column by "_" and extract the first part (The contig name)
  split_values_contig <- sapply(strsplit(merged_data[, 1], "_"), function(x) x[1])
  
  # Split the first column by "_" and extract the second part (position in the contig)
  split_values_pos <- sapply(strsplit(merged_data[, 1], "_"), function(x) x[2])
  
  #Add the the contig and pos columns to the total dataset
  merged_data = cbind(merged_data, split_values_contig, split_values_pos)

  
  # Get unique values
  unique_values <- unique(split_values_contig)
  
  # Write the unique values to the output file without a header
  write.table(unique_values, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  # Write the snp, contig, pos, and LG information
  colnames(merged_data) = c('SNP', 'Inversion', 'LG', 'Contig', 'Pos')
  write.table(merged_data, file = paste0('Full_',output_file), row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  cat(paste("Unique values written to", output_file, "\n"))
}

# Example usage:
# Replace "your_directory" with the actual path to your directory
# Replace "output.txt" with the desired output file name
directory_path <- "Datasets/Diag SNPs parallelogram manu"
output_file_path <- "ParallelogramDiagnosticSNPContigs.txt"

# Call the function
process_text_files(directory_path, output_file_path)
