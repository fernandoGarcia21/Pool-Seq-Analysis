setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')
library(R.utils)
library(dplyr)

#' Translate SNP coordinates from contig to reference genome based on SAM alignment, CIGAR string, and strand.
#'
#' @param contig_pos The position of the SNP in the contig (1-based).
#' @param rname The reference sequence name (RNAME in SAM).
#' @param pos The 1-based leftmost mapping position of the first matching base (POS in SAM).
#' @param cigar The CIGAR string from the SAM file.
#' @param strand "+" for forward strand, "-" for reverse strand.
#' @param contig_length The length of the contig. Required for reverse strand calculations.
#' @return The corresponding position in the reference genome, or NA if the position cannot be translated.
#'
#' @examples
#' translate_snp_coord_strand(10, "chr1", 100, "9M1D5M", "+", 20) # Example with deletion, forward strand
#' translate_snp_coord_strand(10, "chr1", 100, "10M", "+", 10) # Example with perfect match, forward strand.
#' translate_snp_coord_strand(15, "chr1", 100, "10M5I", "+", 20) #Example with insertion, forward strand
#' translate_snp_coord_strand(5, "chr1", 100, "5S10M", "+", 20) # Example with soft clipping, forward strand.
#' translate_snp_coord_strand(10, "chr1", 100, "9M1D5M", "-", 20) # Example with deletion, reverse strand.
#' translate_snp_coord_strand(10, "chr1", 100, "10M", "-", 10) # Example with perfect match, reverse strand.
#' translate_snp_coord_strand(15, "chr1", 100, "10M5I", "-", 20) #Example with insertion, reverse strand
#' translate_snp_coord_strand(5, "chr1", 100, "5S10M", "-", 20) # Example with soft clipping, reverse strand.
translate_snp_coord_strand <- function(contig_pos, rname, pos, cigar, strand, contig_length, p_SNP_name) {
  # Check if inputs are valid
  if (!is.numeric(contig_pos) || contig_pos < 1) {
    stop("contig_pos must be a positive integer.")
  }
  if (!is.numeric(pos) || pos < 1) {
    stop("pos must be a positive integer.")
  }
  if (!is.character(cigar)) {
    stop("cigar must be a character string.")
  }
  if (!strand %in% c("+", "-")) {
    stop("strand must be '+' or '-'.")
  }
  if (!is.numeric(contig_length) || contig_length < 1){
    stop("contig_length must be a positive integer")
  }
  
  cigar_ops <- unlist(regmatches(cigar, gregexpr("[MIDNSHP=X]", cigar)))
  cigar_lens <- as.numeric(unlist(regmatches(cigar, gregexpr("[0-9]+", cigar))))
  
  ref_pos <- pos
  contig_consumed <- 0
  ref_consumed <- 0
  
  if (strand == "-") {
    contig_pos <- contig_length - contig_pos + 1 # Convert to reverse complement coordinates.
    cigar_ops <- rev(cigar_ops)
    cigar_lens <- rev(cigar_lens)
  }
  
  df_extra_info[df_extra_info$SNP == p_SNP_name, ]$strand <<- strand
  for (i in 1:length(cigar_ops)) {
    op <- cigar_ops[i]
    len <- cigar_lens[i]
    if (op == "M" || op == "=" || op == "X") {
      # Match or mismatch
      if (contig_consumed + len >= contig_pos) {
        df_extra_info[df_extra_info$SNP == p_SNP_name, ]$cigar_region <<- op
        if (strand == "+"){
          return(ref_pos + (contig_pos - contig_consumed) - 1)
        } else {
          return(ref_pos + (contig_consumed + len - contig_pos)) #Corrected reverse calculation
        }
      }
      ref_consumed <- ref_consumed + len
      contig_consumed <- contig_consumed + len
      ref_pos <- ref_pos + len
    } else if (op == "I") {
      # Insertion
      if (contig_consumed + len >= contig_pos) {
        df_extra_info[df_extra_info$SNP == p_SNP_name, ]$cigar_region <<- op
        return(NA) # SNP is within an insertion.
      }
      contig_consumed <- contig_consumed + len
    } else if (op == "D" || op == "N") {
      # Deletion or skipped region
      if (contig_consumed + 1 > contig_pos) {
        df_extra_info[df_extra_info$SNP == p_SNP_name, ]$cigar_region <<- op
        if (strand == "+"){
          return(ref_pos + (contig_pos - contig_consumed) - 1)
        } else {
          return(ref_pos + (contig_consumed + len - contig_pos)) #Corrected reverse calculation
        }
      }
      ref_consumed <- ref_consumed + len
      ref_pos <- ref_pos + len
      
    } else if (op == "S") {
      # Soft clipping
      if (contig_consumed + len >= contig_pos) {
        df_extra_info[df_extra_info$SNP == p_SNP_name, ]$cigar_region <<- op
        if (strand == "+"){
          return(pos -1 + (contig_pos - contig_consumed))
        } else {
          return(pos + (contig_consumed + len - contig_pos)) #Corrected reverse calculation
        }
      }
      contig_consumed <- contig_consumed + len
      
    } else if (op == "H" || op == "P") {
      # Hard clipping or padding
      # These operations do not consume contig or reference bases.
    } else {
      stop(paste("Unknown CIGAR operation:", op))
    }
  }
  
  return(NA) # SNP position exceeds contig length or is otherwise unmappable.
}


# Function to check if a flag bit is set
is_flag_set <- function(flag, bit) {
  
  # Get the bit (using bitwAnd())
  # Bit position (0-based, rightmost is 0) so we subtract 1
  #bit_value = bitwShiftR(flag, bit-1) %>% bitwAnd(1)
  
  return((flag %>% bitwAnd(bit)) != 0)
}



file_old_inv_diagnostic_snps = 'Figures/Inversions/Full_ParallelogramDiagnosticSNPContigs.trimmed.200.txt';
file_sam_inversions = 'Figures/Inversions/Parallelogram.Inversions.trimmed.200.mappedContigs.cigar.sam';
file_contig_lengths = 'Figures/Inversions/InversionsParallelogram.contigLengths.trimmed.200.txt';

data_old_inv_diagnostic_snps = read.table(file_old_inv_diagnostic_snps, header = TRUE, check.names = FALSE)
data_contig_lengths = read.table(file_contig_lengths, header = FALSE, check.names = FALSE)
colnames(data_contig_lengths) = c('Contig','Length')

#Load only specific fields of the sam file, for better performance, you can do this.
#sam_fields = c("qname", "flag", "rname", "pos", "cigar", "mapq")
#param = ScanBamParam(what = sam_fields)
#sam_partial = scanBam(file_sam_inversions, param = param)


sam_lines = readLines(file_sam_inversions)
filtered_lines = c()


#Create a DF from the SAM file for convenience
df_sam = data.frame()

for (line in sam_lines) {
  if (startsWith(line, "@")) {
    filtered_lines <- c(filtered_lines, line)
    next
  }

  fields <- strsplit(line, "\t")[[1]]
  df_sam = rbind(df_sam, data.frame(qname = fields[1],
                                    rname = fields[3],
                                    cigar = fields[6], 
                                    quality = fields[5],
                                    flag = as.integer(fields[2]), 
                                    pos = as.integer(fields[4]) - 1 # 1-based leftmost mapping POSition
                                    ))

 
}

#Exclude alignments with quality smaller than 20
df_sam = df_sam[df_sam$quality >= 20,]

#Keep only alignments where LG info exists (i.e. exclude scaffolds)
df_sam = df_sam[grepl("^LG", df_sam$rname), ]
df_sam$LG = sapply(strsplit(df_sam$rname, "_"), function(x) x[1])

#Exclude alignments where the LG info in the reference field has not the format 'LG12_SUPER_3'
# Count the number of underscores in each rname
underscore_counts = sapply(strsplit(df_sam[['rname']], "_"), length) - 1

# Filter the data frame to keep rows with two underscores
df_sam = df_sam[underscore_counts == 2, ]

#############################################################################################
#Add 'LG' to the number of the LG
#data_old_inv_diagnostic_snps$LG = paste0('LG',data_old_inv_diagnostic_snps$LG)
#############################################################################################

#Store the translated coordinates, then append this array to the SNPs DF
translated_coordinates_array = c()
translated_LG_Super_array = c()
df_extra_info = data.frame(SNP = data_old_inv_diagnostic_snps[,'SNP'], strand = NA, cigar_region = NA)

for(i_snp in 1:dim(data_old_inv_diagnostic_snps)[1]){
  tmp_translated_coordinate = NA
  tmp_translated_LG_Super = NA
  tmp_snp = data_old_inv_diagnostic_snps[i_snp,]
  tmp_contig_length = data_contig_lengths[data_contig_lengths$Contig == tmp_snp$Contig, ]$Length
  
  #Exctract the alignment information of the SNP
  tmp_sam_snp = df_sam[df_sam$qname == tmp_snp$SNP & df_sam$LG == tmp_snp$LG, ]
  
  if(dim(tmp_sam_snp)[1] > 0){
    #Find the main alignment based on the flag and discard the supplementary alignments
    #The primary alignment has the flag 0x100 (256) bit unset
    
    tmp_sam_snp = tmp_sam_snp %>%
      filter(!is_flag_set(flag, 4)) %>% # Exclude unmapped
      filter(!is_flag_set(flag, 256)) %>% # Exclude secondary
      filter(!is_flag_set(flag, 2048)) # Exclude duplicates
    
    
    #print(paste('There are ', dim(tmp_sam_snp)[1], 'Good Alignments'))
    
    #Find the alignment with the largest number of M
    if(dim(tmp_sam_snp)[1] == 1){
      
      tmp_index_main_alignment = 1
      #Use the main alignment to translate the SNP coordinate
      tmp_chosen_alignment = tmp_sam_snp[tmp_index_main_alignment,]
      tmp_translated_LG_Super = tmp_chosen_alignment$rname
      cigar = tmp_chosen_alignment$cigar
      ref_start = tmp_chosen_alignment$pos
      flag = tmp_chosen_alignment$flag
      strand = ''
      
      #The bitwise AND operation flag & 16 checks if the 16th bit of the flag is set.
      #if the 16th bit is not set (i.e., the result is 0), it indicates a forward strand alignment.
      if (!is_flag_set(flag, 16)) { # Forward strand
        strand = '+'
        #print(paste('Contig', tmp_snp$Contig, 'Index ',i_snp ))
        #print(cigar)
        
      } else { # Reverse strand
        strand = '-'
      }
      
      #Call the function to translate the coordinate of the SNP
      tmp_translated_coordinate = as.numeric(translate_snp_coord_strand(contig_pos = tmp_snp$Pos, tmp_chosen_alignment$rname, ref_start, cigar, strand, tmp_contig_length, tmp_snp$SNP))
      
      if(!is.na(tmp_translated_coordinate)){
        #The cigar positions are 0-based, so we must add 1 to get get the SNP position of the Pool-Seq SNPs
        tmp_translated_coordinate = tmp_translated_coordinate +1
      }
      
    }else{
      print(paste('The contig', tmp_snp$Contig, 'has multiple good alignments, so excluded.'))
    }
    
  }#If alignments were found for the SNP
  else{
    #print(paste('The contig', tmp_snp$Contig, 'has no alignments in the SAM file after the quality filters.'))
  }
  
  #Append the translated coordinate to the array
  translated_coordinates_array = append(translated_coordinates_array, tmp_translated_coordinate)
  translated_LG_Super_array = append(translated_LG_Super_array, tmp_translated_LG_Super)
  
}


data_old_inv_diagnostic_snps = cbind(data_old_inv_diagnostic_snps, Translated_pos = translated_coordinates_array, LG_SUPER = translated_LG_Super_array, strand = df_extra_info$strand, cigar_region = df_extra_info$cigar_region)
print(paste(dim(data_old_inv_diagnostic_snps[!is.na(data_old_inv_diagnostic_snps$Translated_pos), ])[1], 'SNPs sucessfully translated'))

write.table(file='CorrectedTrimmed.TranslatedSNPsFromParallelogram.200.txt', data_old_inv_diagnostic_snps, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')



