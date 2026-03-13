# Make the gene_to_go file for topGO

# load in file containing ProteinIDs and corresponding GO terms
df <- read.csv("~/Desktop/Stomentosus_proteinID_GOterm.csv", header = T)

# functions to manipulate df to show porteinID and all the GO terms assigned to it
merge_rows <- function(data, identifier_col, merge_col) {
  # Group the data by the identifier column
  grouped_data <- aggregate(data[[merge_col]], by = list(data[[identifier_col]]), FUN = paste, collapse = ",")
  
  # Rename the columns of the grouped data
  colnames(grouped_data) <- c(identifier_col, merge_col)
  
  # Return the merged data frame
  return(grouped_data)
}

merged_data <- merge_rows(df, "ProteinID", "GOTerm")
print(merged_data)

# this is the GO mapping file - write it
write.table(merged_data, file = "~/Desktop/Stomentosus_gene_to_go.txt", sep = "\t",
            row.names = FALSE)
