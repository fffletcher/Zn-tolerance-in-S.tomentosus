# fungal DE gene Expression - BIOASSAY
# fungal inocualtaion at 0 mM - ZnT vs ZnS

# load libraries
library(tidyverse)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)
library(ggplot2)
library(scales)
library(ggvenn)

ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Transcriptome Analysis/counts.txt",
                                    header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# make a df for tissue from 220
counts <- ref230_ %>%
  dplyr::select(grep("bioassay_0Zn", colnames(ref230_), value = T))
colnames(counts)

# exclude 230 rep 4 - outlier
counts <- counts %>%
  dplyr::select(!grep("X230_bioassay_0Zn_Rep4", colnames(counts), value = T)) 

# make col data files from samples in the df
sample_names <- colnames(counts)

# Function to extract isolate from the name
extract_info <- function(sample) {
  # Extract treatment (digits before 'Zn')
  treatment <- sub(".*_([0-9]+\\.?[0-9]*)Zn_.*", "\\1", sample)
  isolate <- ifelse(grepl("control_pine", sample), "none", sub(".*X(\\d+).*", "\\1", sample))
  inoculated <- ifelse(grepl("control_pine", sample), "none", "inoculated")
  return(c(sample, treatment, isolate, inoculated))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("sample", "treatment", "isolate", "inocualted")
row.names(col_data) <- NULL
col_data$treatment <- as.factor(col_data$treatment)
print(col_data)

#remove any NA datasets (change to 0s)
counts[is.na(counts)] <- 0

# Create DESeq object (note design only by treatment)
ddsMat <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = col_data,
                                 design = ~ isolate)

# Filtering: removing rows of the DESeqDataSet that have no counts, 
# or only a single count across all samples.
keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

# rlog transformation
# variance stabilizing transformation (VST)
vst <- vst(ddsMat, blind = FALSE)

# generate a PCA based of the DESeq object
plotPCA(vst, intgroup = c("isolate")) +
  aes(shape = isolate, color = isolate) #+ geom_text(label = colnames(vst))

# Run the DESeq stats - gives pairwise comparisons
# this compares all treatments to the 0 baseline and stores in 1 df
dds <- DESeq(ddsMat)
results_names <- resultsNames(dds)
results_names <- results_names[-1] # removes the "Intercept" col - not needed
results_names

# make a for loop to make a df of each comparison

# create a template df to fill in with the for loop
# df for logFC results.
res_df <- data.frame("gene_id" = character(), 
                     "baseMean" = numeric(), 
                     "log2FoldChange" = numeric(), 
                     "lfcSE" = numeric(),
                     "stat" = numeric(),
                     "pvalue" = numeric(),
                     "padj" = numeric(),
                     "gene_name" = character(),
                     "result_name" = character())

# also create a "shrunken" df 
# logFC this normalizes low expressed genes to be less significant

res_shrunken_df <- data.frame("gene_id" = character(), 
                              "baseMean" = numeric(), 
                              "log2FoldChange" = numeric(), 
                              "lfcSE" = numeric(),
                              "stat" = numeric(),
                              "pvalue" = numeric(),
                              "padj" = numeric(),
                              "gene_name" = character(),
                              "result_name" = character())

# for loop to go through the list of results_names
# and populate both dfs above with the results from each comparison

for(i in 1:length(results_names)) {
  # grabbing the name of the result file i
  results_name <- results_names[i]
  # populating the res_df with results(dds)
  res <- results(dds, name = results_name)
  # populating res shrunken lfc with lfcShrink
  res_shrunken <- lfcShrink(dds, coef = results_name,  res = res)
  # populating data.frame 1 : temp_res_df
  tmp_res_df <- res %>% as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    mutate(result_name = results_name)
  
  # populating data.frame 1 : temp_res_shrunken
  tmp_res_shrunken_df <- res_shrunken %>% as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    mutate(result_name = results_name)
  
  # Append to full data.frame
  res_df <- bind_rows(res_df, tmp_res_df)
  res_shrunken_df <- bind_rows(res_shrunken_df, tmp_res_shrunken_df)
}

# use shrunken
res_df
res_shrunken_df

# write these to csv and save for further anaysis
#write.csv(res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/pine_0_inocualted_v_not.csv")

# Let's make list of all SIG gnees with  P <0.05 and FC > 1.5 in the comparisons
sig_genes <- res_shrunken_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.5)

sum(sig_genes$log2FoldChange >= 1.5)
sum(sig_genes$log2FoldChange <= -1.5)

# make a volcano plot of the differentiall expressed genes
keyvals <- ifelse(
  res_shrunken_df$log2FoldChange <= -1.5 & res_shrunken_df$padj < 0.05, '#E69F00',
  ifelse(res_shrunken_df$log2FoldChange >= 1.5 & res_shrunken_df$padj < 0.05, '#0072B2',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#0072B2'] <- 'upregulated'
names(keyvals)[keyvals == '#E69F00'] <- 'down-reguated'

EnhancedVolcano(res_shrunken_df,
                lab = rownames(res_shrunken_df),
                selectLab = c(''),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-20, 20))


# GO enrichment
library(topGO)

# read the GO mapping file - with IDsep specification
gene_to_go_mapping <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Stomentosus_gene_to_go.txt", 
                                   sep = "\t", IDsep = ",")

# 230 tissue assay RNA comparing 0 and 10 mM Zn
# the naming conventions of our gene_id column are strange
# we want only the protein ID in the gene-id column

# Extract the protein ID and update the gene_id column
# uses dfs from the first DEG script 
de_genes <- res_shrunken_df %>%
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
#write.csv(sig_go_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_DEG_GO_10mM.csv")

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