# find go terms in a Go Term
# this script uses the GeneList that was generated from p value and l2fc
# will need to edit if want to look at only genes significant (not l2fc threshold)

# Replace 'GO:XXXXXXX' with your actual GO term ID
go_term_id <- "GO:0016788"

# Get all genes annotated to this GO term
genes_for_term <- genesInTerm(my_go_dataMF, whichGO = go_term_id)[[1]]

# Get the genes that were selected (i.e., marked as 1 in geneList)
selected_genes <- names(geneList)[as.numeric(as.character(geneList)) == 1]

# Intersect with genes annotated to the GO term
contributing_genes <- intersect(genes_for_term, selected_genes)

contributing_genes

