library(topGO)
library(dplyr)
library(ggplot2)

# read the GO mapping file - with IDsep specification
gene_to_go_mapping <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Stomentosus_gene_to_go.txt", 
                                   sep = "\t", IDsep = ",")

# 230 tissue assay RNA comparing 0 and 10 mM Zn
# the naming conventions of our gene_id column are strange
# we want only the protein ID in the gene-id column

# Extract the protein ID and update the gene_id column
# uses dfs from the first DEG script 
de_genes <- cont_v_10_res_shrunken_df %>%
  mutate(gene_id = sub(".*proteinId=([0-9]+).*", "\\1", gene_id))

# all DE genes in the experimnent = Universe


# first, get only the name and the p value
deGenes_Universe <- de_genes %>%
  dplyr::select(gene_id, padj)

# sort by the membership value (in descending order)
deGenes_Universe <- deGenes_Universe[order(-deGenes_Universe$padj), ]

# convert to topGO's genelist format
topgoUniverse <- deGenes_Universe$padj
names(topgoUniverse) <- deGenes_Universe$gene_id
topgoUniverse <- na.omit(topgoUniverse)

# Gene Selection needs to run off a function 
# This function returns the genes with p value under 0.05
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)
}

# This is where we create the topGO object
# This needs to be repeated with each GO Term category
# Molecular Function  = "MF"
# Cellular Component = "CC"
# Biological Pathway = "BP"
my_go_dataBP <- new("topGOdata",
                    ontology    = "BP",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

## Some stats to look at which GO terms are significantly over represented
## in the significantly differentially expressed genes

## Calculate fisher exact test using 'weight01' algorithm:
resultTableBP <- runTest(my_go_dataBP, algorithm = "weight01", statistic = "fisher")
allResBP <- GenTable(my_go_dataBP, 
                     weightFisher = resultTableBP, 
                     topNodes = 20)
allResBP

# CC
my_go_dataCC <- new("topGOdata",
                    ontology    = "CC",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

## Some stats to look at which GO terms are significantly over represented
## in the significantly differentially expressed genes

## Calculate fisher exact test using 'weight01' algorithm:
resultTableCC <- runTest(my_go_dataCC, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(my_go_dataCC, 
                     weightFisher = resultTableCC, 
                     topNodes = 20)
allResCC

# MF
my_go_dataMF <- new("topGOdata",
                    ontology    = "MF",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

## Some stats to look at which GO terms are significantly over represented
## in the significantly differentially expressed genes

## Calculate fisher exact test using 'weight01' algorithm:
resultTableMF <- runTest(my_go_dataMF, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(my_go_dataMF, 
                     weightFisher = resultTableMF, 
                     topNodes = 20)
allResMF

# plot bar plot
# Load required libraries

# Add category labels
allResBP$Category <- "Biological Process"
allResMF$Category <- "Molecular Function"
allResCC$Category <- "Cellular Component"

# Combine all data
go_df <- bind_rows(allResBP, allResMF, allResCC)

# Filter for significant terms (adjusted p-value < 0.05)
sig_go_df <- go_df %>% filter(weightFisher < 0.05)
write.csv(sig_go_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_DEG_GO_10mM.csv")

# Optional: reorder terms by p.adjust for better plotting
sig_go_df <- sig_go_df %>%
  mutate(Term = fct_reorder(Term, weightFisher))
sig_go_df$weightFisher <- as.numeric(sig_go_df$weightFisher)


# Combine Term and GO.ID for y-axis labels
sig_go_df <- sig_go_df %>%
  mutate(Label = paste0(Term, " (", GO.ID, ")")) %>%
  mutate(Label = fct_reorder(Label, Significant))

# Plot
ggplot(sig_go_df, aes(x = Significant, y = Label, fill = weightFisher)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0(Significant, " / ", round(Expected, 1))), hjust = -0.1, size = 3) +
  facet_grid(Category ~ ., scales = "free", space = "free") +
  scale_fill_gradient(low = "skyblue", high = "darkblue", name = "Fisher p-value") +
  labs(
    title = "Significant GO Terms by Category",
    x = "Number of Significant Genes",
    y = "GO Term (GO:ID)"
  ) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text = element_text(size = 10),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey90", color = NA)) +
  xlim(0,15)
                            
# pull genes contributing to significant go terms and map expression over Zn gradient
# uses res_shrunken_df from the first .R script

# Get genes associated with the GO term from topGO
go_term <- "GO:0043168"
# here change data to the crrroect topGO object = my_go_dataMF, my_go_dataBP or my_go_dataCC
go_genes <- genesInTerm(my_go_dataMF, go_term)[[1]]  

deg_data <- res_shrunken_df
deg_data$gene_id <- sub(".*proteinId=([0-9]+).*", "\\1", deg_data[[1]])
deg_data$result_name <- factor(deg_data$result_name, 
                               levels = c("treatment_0.1_vs_0",
                                          "treatment_1_vs_0", 
                                          "treatment_2.5_vs_0",
                                          "treatment_10_vs_0"))

# Filter DEG data for these genes
deg_goi <- deg_data %>%
  filter(gene_id %in% go_genes) %>%
  mutate(significance = ifelse(padj < 0.05, "Significant", "Not Significant"))

# Plot expression across treatments
ggplot(deg_goi, aes(x = result_name, y = log2FoldChange, group = gene_id, color = gene_id)) +
  geom_line() +
  geom_point(aes(shape = significance, alpha = significance)) +
  scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.4)) +
  labs(title = paste("Expression of Genes in", go_term),
       x = "Treatment",
       y = "Log2 Fold Change vs Control") +
  theme_minimal() +
  theme(legend.position = "none")



