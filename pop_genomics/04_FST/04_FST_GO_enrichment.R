# Load required libraries
library(topGO)
library(dplyr)
library(ggplot2)
library(tidyr)

# read in your gene-to-GO mapping
geneID2GO <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Stomentosus_gene_to_go.txt", 
                          sep = "\t", IDsep = ",")

# read in list of genes from FST analysis
# use this one for genes hit by SNPs in top 1% FST windows
# fst_snps <- read.csv("~/Desktop/Stomentosus_PopGen/FST/FST_top_1p_annotated_snps_with_genes.csv", header = TRUE)
# fst_snps <- fst_snps[!fst_snps$WEIR_AND_COCKERHAM_FST == 0,] # remove anything w/ FST = 0
# fst_snps$Gene <- as.character(fst_snps$Gene) # for later

# use this one for top 5%
fst_snps <- read.csv("~/Desktop/Stomentosus_PopGen/FST/FST_top_5p_annotated_snps_with_genes.csv", header = TRUE)
fst_snps <- fst_snps[!fst_snps$WEIR_AND_COCKERHAM_FST == 0,] # remove anything w/ FST = 0
fst_snps$Gene <- as.character(fst_snps$Gene) # for later

# use this one for genes hit by SNPs in only FST = 1 windows
#genesOfInterest <- read.csv("~/Desktop/Stomentosus_PopGen/FST/FST1_annotated_snps_with_genes.csv", header = TRUE)

# count SNPs per gene and region
snp_counts <- fst_snps %>%
  group_by(Gene, Annotation) %>%
  summarise(SNP_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Annotation, values_from = SNP_count, values_fill = 0)

# add total SNPs per gene
snp_counts <- snp_counts %>%
  mutate(Total_SNPs = rowSums(across(c(Promotor, Genic), ~replace_na(., 0)))) 

# Extract the gene ID column (protein ID in this case)
genesOfInterest <- fst_snps$Gene
genesOfInterest <- unique(genesOfInterest)
genesOfInterest <- na.omit(genesOfInterest)

# how many genes in list - plus how many are unannotated
total_genes_in_list <- length(genesOfInterest)
annotated_genes <- genesOfInterest[genesOfInterest %in% names(geneID2GO)]
num_annotated <- length(annotated_genes)
num_unannotated <- total_genes_in_list - num_annotated
cat("Total genes in FST list:", total_genes_in_list, "\n")
cat("Annotated genes (with GO terms):", num_annotated, "\n")
cat("Unannotated genes (no GO terms):", num_unannotated, "\n")

# all genes are all the genes with GO annotations
allGenes <- names(geneID2GO)
# geneList is a vector with the gene IDs as names
# and a 1/0 if genes are present in the FST list
geneList <- factor(as.integer(allGenes %in% genesOfInterest))
names(geneList) <- allGenes

# create topGO object - Molecular Function
MFdata <- new("topGOdata",
              ontology = "MF",  
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

# run enrichment test (Fisher's exact test)
resultFisherMF <- runTest(MFdata, algorithm = "weight01", statistic = "fisher")
# put in table
resultsTableMF <- GenTable(MFdata,
                         classicFisher = resultFisherMF,
                         orderBy = "classicFisher",
                         topNodes = 10)  # number of terms returned in table
resultsTableMF

# Add contributing genes to results table
resultsTableMF$contributing_genes <- sapply(resultsTableMF$GO.ID, function(go_id) {
  term_genes <- genesInTerm(MFdata, go_id)[[1]]
  hits <- term_genes[term_genes %in% genesOfInterest]
  paste(hits, collapse = ", ")
})

resultsTableMF

# filter for significant GO terms (p < 0.05)
significant_termsMF <- resultsTableMF %>%
  filter(as.numeric(classicFisher) < 0.05)

# Ssplit contributing genes into individual rows
gene_go_pairsMF <- significant_termsMF %>%
  dplyr::select(GO.ID, Term, contributing_genes) %>%
  separate_rows(contributing_genes, sep = ",\\s*")

# join SNP annotation data
detailed_tableMF <- gene_go_pairsMF %>%
  left_join(fst_snps, by = c("contributing_genes" = "Gene")) %>%
  dplyr::select(GO.ID, Term, contributing_genes, SNP, Annotation, WEIR_AND_COCKERHAM_FST)

detailed_tableMF

summary_tableMF <- detailed_tableMF %>%
  group_by(GO.ID, Term, contributing_genes, Annotation) %>%
  summarise(SNP_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Annotation, values_from = SNP_count, values_fill = 0) %>%
  mutate(Total_SNPs = rowSums(across(c("Promotor", "Genic"), ~replace_na(., 0))))

summary_tableMF

# create topGOdata object - Biological Process
BPdata <- new("topGOdata",
              ontology = "BP",  
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

resultFisherBP <- runTest(BPdata, algorithm = "weight01", statistic = "fisher")

resultsTableBP <- GenTable(BPdata,
                         classicFisher = resultFisherBP,
                         orderBy = "classicFisher",
                         topNodes = 10) 
resultsTableBP

resultsTableBP$contributing_genes <- sapply(resultsTableBP$GO.ID, function(go_id) {
  term_genes <- genesInTerm(BPdata, go_id)[[1]]
  hits <- term_genes[term_genes %in% genesOfInterest]
  paste(hits, collapse = ", ")
})

resultsTableBP

significant_termsBP <- resultsTableBP %>%
  filter(as.numeric(classicFisher) < 0.05)

gene_go_pairsBP <- significant_termsBP %>%
  dplyr::select(GO.ID, Term, contributing_genes) %>%
  separate_rows(contributing_genes, sep = ",\\s*")

detailed_tableBP <- gene_go_pairsBP %>%
  left_join(fst_snps, by = c("contributing_genes" = "Gene")) %>%
  dplyr::select(GO.ID, Term, contributing_genes, SNP, Annotation, WEIR_AND_COCKERHAM_FST)

detailed_tableBP

summary_tableBP <- detailed_tableBP %>%
  group_by(GO.ID, Term, contributing_genes, Annotation) %>%
  summarise(SNP_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Annotation, values_from = SNP_count, values_fill = 0) %>%
  mutate(Total_SNPs = rowSums(across(c("Promotor", "Genic"), ~replace_na(., 0))))

summary_tableBP

# create topGOdata object - Cellular Component
CCdata <- new("topGOdata",
              ontology = "CC",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)

resultFisherCC <- runTest(CCdata, algorithm = "weight01", statistic = "fisher")

resultsTableCC <- GenTable(CCdata,
                         classicFisher = resultFisherCC,
                         orderBy = "classicFisher",
                         topNodes = 10) 
resultsTableCC

resultsTableCC$contributing_genes <- sapply(resultsTableCC$GO.ID, function(go_id) {
  term_genes <- genesInTerm(CCdata, go_id)[[1]]
  hits <- term_genes[term_genes %in% genesOfInterest]
  paste(hits, collapse = ", ")
})

significant_termsCC <- resultsTableCC %>%
  filter(as.numeric(classicFisher) < 0.05)

gene_go_pairsCC <- significant_termsCC %>%
  dplyr::select(GO.ID, Term, contributing_genes) %>%
  separate_rows(contributing_genes, sep = ",\\s*")

detailed_tableCC <- gene_go_pairsCC %>%
  left_join(fst_snps, by = c("contributing_genes" = "Gene")) %>%
  dplyr::select(GO.ID, Term, contributing_genes, SNP, Annotation, WEIR_AND_COCKERHAM_FST)

detailed_tableCC

summary_tableCC <- detailed_tableCC %>%
  group_by(GO.ID, Term, contributing_genes, Annotation) %>%
  summarise(SNP_count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Annotation, values_from = SNP_count, values_fill = 0) %>%
  mutate(Total_SNPs = rowSums(across(c("Promotor", "Genic"), ~replace_na(., 0))))

summary_tableCC


## combine the results table and output csv
bp <- resultsTableBP[resultsTableBP$classicFisher<0.05,]
cc <- resultsTableCC[resultsTableCC$classicFisher<0.05,]
mf <- resultsTableMF[resultsTableMF$classicFisher<0.05,]

bp$Ontology <- "Biological Process"
cc$Ontology <- "Cellular Component"
mf$Ontology <- "Molecular Function"

table1 <- rbind(bp, cc, mf)
write.csv(table1, "~/Desktop/Stomentosus_PopGen/FST/top5p_GOterms.csv")



## plotting GO data
## rename these for whichever group of terms you're working with

resultsTable <- resultsTableMF
summary_table <- summary_tableMF
detailed_table <- detailed_tableMF


resultsTable <- resultsTable %>%
  mutate(log_p = -log10(as.numeric(classicFisher)),
         significant = as.numeric(classicFisher) < 0.05) %>%
  arrange(classicFisher)

# Plot the top 20 terms and highlight significant ones
ggplot(resultsTable, aes(x = reorder(Term, log_p), y = log_p, fill = significant)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                    labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
                    name = "Significance") +
  labs(title = "Top 10 Enriched GO Terms (Cellular Component)",
       x = "GO Term",
       y = expression(-log10)) +
  theme_minimal()


# plot the contributing genes for each significant go term

# Loop through each GO term
unique_terms <- unique(summary_table$Term)

for (term in unique_terms) {
  # Filter for the current GO term
  term_data <- summary_table %>%
    filter(Term == term) %>%
    dplyr::select(contributing_genes, Genic, Promotor)
  
  # Reshape for stacked bar plot
  snp_long <- term_data %>%
    pivot_longer(cols = c(Genic, Promotor), names_to = "SNP_Type", values_to = "Count")
  
  # Create the plot
  p <- ggplot(snp_long, aes(x = contributing_genes, y = Count, fill = SNP_Type)) +
    geom_bar(stat = "identity") +
    labs(title = paste("SNP Distribution for GO Term:", term),
         x = "Contributing Gene",
         y = "SNP Count") +
    scale_fill_manual(values = c("Genic" = "firebrick", "Promotor" = "steelblue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}


# fst distribution across GO terms

# want this in x axis order of significance 
# create an ordered list of GO terms based on p-value from resultsTable
# this is to set the order of GO terms by significance on x axis
ordered_terms <- resultsTable %>%
  arrange(classicFisher) %>%
  pull(Term)

# apply this order to the "Term" factor in your detailedTable
detailed_table$Term <- factor(detailed_table$Term, levels = ordered_terms)

# boxplot
ggplot(detailed_table, aes(x = Term, y = WEIR_AND_COCKERHAM_FST)) +
  geom_boxplot(fill = "steelblue") +
  labs(title = "FST Distribution per GO Term",
       x = "GO Term",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

# jitter
ggplot(detailed_table, aes(x = Term, y = WEIR_AND_COCKERHAM_FST)) +
  geom_jitter(width = 0.2, alpha = 1, color = "steelblue") +
  labs(title = "FST Values per GO Term",
       x = "GO Term",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

# combined
ggplot(detailed_table, aes(x = Term, y = WEIR_AND_COCKERHAM_FST)) +
  geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = .5, color = "darkblue") +
  labs(title = "FST Values per GO Term",
       x = "GO Term",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)


# FST per gene in GO term? 
ggplot(detailed_table[detailed_table$Term=="zinc ion binding",], aes(x = contributing_genes, y = WEIR_AND_COCKERHAM_FST)) +
  geom_boxplot(aes(middle = mean(WEIR_AND_COCKERHAM_FST)), fill = "steelblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = .5, color = "darkblue") +
  labs(title = "FST Distribution per gene in GO Term: zinc ion binding ",
       x = "Gene",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

# plot mean instead of median for boxplots - zinc ion binding
# this is to set the order of GO terms by mean on x axis
zn <- detailed_table[detailed_table$Term=="zinc ion binding",]

ordered_terms1 <- zn %>%
  group_by(contributing_genes) %>%
  summarise(median_fst = median(WEIR_AND_COCKERHAM_FST, na.rm = TRUE)) %>%
  arrange(desc(median_fst))

zn$contributing_genes <- factor(zn$contributing_genes,
                                levels = ordered_terms1$contributing_genes)

ggplot(zn, aes(x = contributing_genes, y = WEIR_AND_COCKERHAM_FST)) +
  geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = .5, color = "darkblue") +
  labs(title = "FST Distribution per gene in GO Term: zinc ion binding ",
       x = "Gene",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)


ggplot(zn, aes(x = contributing_genes, y = WEIR_AND_COCKERHAM_FST)) +
  geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue" ) +
  geom_jitter(width = 0.2, alpha = .5, color = "darkblue") +
  labs(title = "FST Distribution per gene in GO Term: zinc ion binding",
       x = "Gene",
       y = "FST Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

  