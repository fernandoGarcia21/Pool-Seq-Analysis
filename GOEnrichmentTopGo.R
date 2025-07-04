setwd('C:/Users/dgarciac/Documents/PhD Project/Third Chapter/PoolSeq-Analyses/')

#BiocManager::install("topGO", force=TRUE)

# Load the ggplot2 library
library(ggplot2)
library(stringr)
library(topGO)

population_name = 'RED'

# 1. Read the GAF file
# We need to skip comment lines (starting with '!')
# sep="\t" indicates tab-delimited
# header=FALSE because there's no header row after the comments
# fill=TRUE handles potential uneven rows if any (though GAF is usually strict)
gaf_file_path = "Annotation/GCF_037325665.1-RS_2024_12_gene_ontology.gaf"

# Find the number of comment lines to skip
# Read the file line by line to count comments
con = file(gaf_file_path, "r")
skip_lines = 0
while (TRUE) {
  line = readLines(con, n = 1)
  if (length(line) == 0) break # End of file
  if (startsWith(line, "!")) {
    skip_lines = skip_lines + 1
  } else {
    break # First non-comment line found
  }
}
close(con)

# Now read the data, skipping the determined number of lines
gaf_data = read.delim(gaf_file_path,
                       sep = "\t",
                       header = FALSE,
                       skip = skip_lines,
                       quote = "", # GAF files typically don't use quotes
                       col.names = c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier",
                                     "GO_ID", "DB_Reference", "Evidence_Code", "With_From",
                                     "Aspect", "DB_Object_Name", "DB_Object_Synonym",
                                     "DB_Object_Type", "Taxon", "Date", "Assigned_By",
                                     "Annotation_Extension", "Gene_Product_Form_ID"),
                       fill = TRUE,
                       check.names = FALSE) # Prevent R from sanitizing column names

# Clean up the temporary file (optional)
#unlink(gaf_file_path)

# Display the first few rows to check
head(gaf_data)

# 2. Extract relevant columns
# We need DB_Object_ID (column 2), GO_ID (column 5), and Aspect (column 9)
parsed_annotations = gaf_data[, c("DB_Object_ID", "GO_ID", "Aspect")]

# 3. Create the geneID2GO mapping for topGO
# The 'geneID2GO' object for topGO should be a list where:
# - Each element is a character vector of GO IDs.
# - The names of the list elements are the gene IDs.

# Aggregate GO terms by gene ID
geneID2GO_list = tapply(parsed_annotations$GO_ID,
                         parsed_annotations$DB_Object_ID,
                         function(x) as.character(x))



######################### GO Enrichment Analysis ############################

# 1. Your list of genes of interest (foreground genes)
# Example: Using a small set of dummy gene IDs
con = file(paste0("Annotation/",population_name,".geneIds.txt"), "r")
myInterestingGenes = readLines(con)
close(con)



# 2. Your gene universe (background genes)
# Example: A larger set of dummy gene IDs

gtf_file_path = "Annotation/GCF_037325665.1_US_GU_Lsax_2.0_genomic.gtf"

# Find the number of comment lines to skip
# Read the file line by line to count comments
con = file(gtf_file_path, "r")
skip_lines = 0
while (TRUE) {
  line = readLines(con, n = 1)
  if (length(line) == 0) break # End of file
  if (startsWith(line, "#")) {
    skip_lines = skip_lines + 1
  } else {
    break # First non-comment line found
  }
}
close(con)

# Now read the data, skipping the determined number of lines
gtf_data = read.delim(gtf_file_path,
                      sep = "\t",
                      header = FALSE,
                      skip = skip_lines,
                      quote = "", # GAF files typically don't use quotes
                      col.names = c("seaname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
                      fill = TRUE,
                      check.names = FALSE) # Prevent R from sanitizing column names

# Clean up the temporary file (optional)
#unlink(gtf_file_path)

# Display the first few rows to check
head(gtf_data)

# 2. Extract relevant columns
# We need DB_Object_ID (column 2), GO_ID (column 5), and Aspect (column 9)
gtf_attributes = gtf_data[, c("feature", "attribute")]
gtf_attributes = gtf_attributes[gtf_attributes$feature == "gene", ]
gtf_attributes$GeneID <- str_extract(gtf_attributes$attribute, 'GeneID:([0-9]+)')
gtf_attributes$GeneID <- sub("GeneID:", "", gtf_attributes$GeneID)
allGenes = gtf_attributes$GeneID
rm(gtf_attributes, gtf_data)

# You need to define which genes from your 'allGenes' list are in 'myInterestingGenes'
# This creates a factor where 1 indicates an interesting gene, 0 indicates not.
geneList = factor(as.integer(allGenes %in% myInterestingGenes))
names(geneList) <- allGenes


# Define the ontology you want to analyze
# "BP" for Biological Process
# "MF" for Molecular Function
# "CC" for Cellular Component
myOntology = "BP"

# Create the topGOdata object
# allGenes: Your complete background gene set
# geneList: A factor indicating which genes are "interesting"
# annot: Your gene-to-GO term mapping
# nodeSize: Minimum number of genes in a GO term to be considered.
#           Setting this to a small number (e.g., 10) can filter out very specific terms.
GOdata = new("topGOdata",
              ontology = myOntology,
              allGenes = geneList,
              annot = annFUN.gene2GO, # Use annFUN.gene2GO for custom mapping
              gene2GO = geneID2GO_list,
              nodeSize = 10) # Adjust nodeSize as appropriate for your data

#View(geneID2GO_list)

# Run the Fisher's exact test with the 'elim' algorithm
# This algorithm tries to eliminate false positives by considering the GO hierarchy.
resultFisher_elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher")

# Alternatively, 'weight01' is often recommended
resultFisher_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "Fisher")

# You can also run the 'classic' algorithm for comparison (often has more significant terms)
resultFisher_classic <- runTest(GOdata, algorithm = "classic", statistic = "Fisher")




# Generate a table of results
# orderBy: How to order the results (e.g., "elimFisher" for p-value from elim algorithm)
# ranksOf: Which test results to display p-values for
# topNodes: Number of top GO terms to display
go_results_elim <- GenTable(GOdata,
                            elimFisher = resultFisher_elim,
                            classicFisher = resultFisher_classic, # Include classic for comparison
                            orderBy = "elimFisher",
                            ranksOf = "elimFisher",
                            topNodes = 100)

# Display the results
#View(go_results_elim)

# You can also filter by significance
# For example, terms with elimFisher p-value < 0.05
significant_terms <- go_results_elim[as.numeric(go_results_elim$elimFisher) < 0.05, ]
#View(significant_terms)
#write.table(significant_terms, file = paste0("Annotation/",population_name,".Enrichment.tsv"))

# Understanding the columns:
# GO.ID: The Gene Ontology ID
# Term: The name of the GO term
# Annotated: Total number of genes annotated to this term in your background.
# Significant: Number of genes from your interesting gene list annotated to this term.
# Expected: Expected number of interesting genes for this term by chance.
# elimFisher/classicFisher: P-value from the respective statistical test.


# Generate a table of results, including all tests you ran
go_results_table <- GenTable(GOdata,
                             elimFisher = resultFisher_elim,
                             weight01Fisher = resultFisher_weight01, # Add weight01 results
                             classicFisher = resultFisher_classic,   # Include classic for comparison
                             orderBy = "elimFisher", # You can order by elimFisher or weight01Fisher
                             ranksOf = "elimFisher",
                             topNodes = length(usedGO(GOdata))) # Get all calculated terms for adjustment




# Ensure the p-value column is numeric before adjustment
go_results_table$elimFisher <- as.numeric(go_results_table$elimFisher)
go_results_table$weight01Fisher <- as.numeric(go_results_table$weight01Fisher)
go_results_table$classicFisher <- as.numeric(go_results_table$classicFisher)

# --- APPLY FDR CORRECTION (Benjamini-Hochberg) ---

# For elimFisher p-values
go_results_table$elimFisher_FDR <- p.adjust(go_results_table$elimFisher, method = "BH")

# For weight01Fisher p-values (recommended to use this one as it's often more robust)
go_results_table$weight01Fisher_FDR <- p.adjust(go_results_table$weight01Fisher, method = "BH")

# For classicFisher p-values (for comparison, usually has many significant terms without hierarchy consideration)
go_results_table$classicFisher_FDR <- p.adjust(go_results_table$classicFisher, method = "BH")


# Display the results with FDR adjusted p-values
View(go_results_table)

#Write the FDR test restults
write.table(go_results_table, file = paste0("Annotation/",population_name,".FDR.Enrichment.tsv"))

# You can now filter by the FDR adjusted p-values (q-values)
# For example, terms with elimFisher_FDR < 0.05
significant_terms_FDR_elim <- go_results_table[go_results_table$elimFisher_FDR < 0.05, ]
View(significant_terms_FDR_elim)

# Or, if using weight01 (often preferred for its balance)
significant_terms_FDR_weight01 <- go_results_table[go_results_table$weight01Fisher_FDR < 0.05, ]
View(significant_terms_FDR_weight01)


# Save the full table with all p-values and FDR adjusted values
write.table(go_results_table,
            file = paste0("Annotation/",population_name,".Enrichment_with_FDR.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# You might want to save only the significant ones as well
write.table(significant_terms_FDR_weight01, # Or significant_terms_FDR_elim
            file = paste0("Annotation/",population_name,".Significant_Enrichment_FDR.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)






# --- Correct way to find the genes that were actually used (feasible genes) ---
feasible_gene_names <- genes(GOdata)

# Find genes from your initial 'allGenes' that are NOT in 'feasible_gene_names'
genes_not_in_ontology <- setdiff(myInterestingGenes, feasible_gene_names)

# Display the count and the list of unannotated genes
cat("Total genes in your universe (allGenes):", length(allGenes), "\n")
cat("Feasible genes (genes with GO annotations in", myOntology, "):", length(feasible_gene_names), "\n")
cat("Genes not found in the ontology (or below nodeSize):", length(genes_not_in_ontology), "\n")
cat("Genes found in the ontology from my interest:", length(myInterestingGenes) - length(genes_not_in_ontology), "\n")




###################################################################
#Check how many genes are annotated for any term
# Get all unique gene IDs that have at least one GO annotation in your GAF file
# (This includes genes annotated to BP, MF, or CC terms)
genes_with_any_go_annotation <- names(geneID2GO_list)

# Now, compare this to your full universe of genes
# (your 'allGenes' object that you derived from the GTF file)

# Count total genes in your universe
total_genes_in_universe <- length(allGenes)

# Count genes with any GO annotation
count_genes_with_any_go <- length(genes_with_any_go_annotation)

# Find genes from your universe that have NO GO annotation at all (across any ontology)
genes_completely_unannotated <- setdiff(allGenes, genes_with_any_go_annotation)
count_completely_unannotated <- length(genes_completely_unannotated)


cat("\n--- Annotation Summary (Any Ontology) ---\n")
cat("Total genes in your universe (from GTF):", total_genes_in_universe, "\n")
cat("Genes with at least ONE GO annotation (BP, MF, or CC):", count_genes_with_any_go, "\n")
cat("Genes with NO GO annotation whatsoever:", count_completely_unannotated, "\n")
cat("Percentage of genes unannotated:",
    round((count_completely_unannotated / total_genes_in_universe) * 100, 2), "%\n")

