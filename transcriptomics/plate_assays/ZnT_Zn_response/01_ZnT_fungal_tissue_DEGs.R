# Stomentosus DE gene Expression
# using reads aligned to 230 genome (better genome)
# What comparisons are wanted? 

# load libraries
library(tidyverse)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)


# set seed
set.seed(1916)

# First look at tissue - look at 230
# counts sheet

ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# make a df for tissue from 230
counts230all <- ref230_ %>%
  dplyr::select(grep("X230_tissue", colnames(ref230_), value = T))
colnames(counts230all)


# remove outliers - identified previously
tissue230counts <- counts230all %>%
  dplyr::select(!grep("_0Zn_Rep4", colnames(counts230all), value = T)) %>%
  dplyr::select(!grep("_1Zn_Rep3", colnames(counts230all), value = T))
# remove all 5 mM conditions
#tissue230counts <- tissue230counts %>%
  #dplyr::select(!grep("_5Zn_Rep", colnames(tissue230counts), value = T))
# results are slightly changed when we do not include 5 mM in the analysis
# if we exclude - we get more genes significant in other treatments!!!
# leave it in for now
colnames(tissue230counts)

# make col data files from samples in the df
sample_names <- colnames(tissue230counts)

# Function to extract treatment from the name
extract_info <- function(sample) {
  # Extract treatment (digits before 'Zn')
  treatment <- sub(".*_([0-9]+\\.?[0-9]*)Zn_.*", "\\1", sample)
  
  return(c(sample, treatment))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("", "treatment")
row.names(col_data) <- NULL
col_data$treatment <- as.factor(col_data$treatment)
print(col_data)

#remove any NA datasets (change to 0s)
tissue230counts[is.na(tissue230counts)] <- 0
# Create DESeq object (note design only by treatment as they are all the same isolate)
ddsMat <- DESeqDataSetFromMatrix(countData = tissue230counts,
                                 colData = col_data,
                                 design = ~ treatment)

# Filtering: removing rows of the DESeqDataSet that have no counts, 
# or only a single count across all samples.
keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

# rlog transformation
# variance stabilizing transformation (VST)
vst <- vst(ddsMat, blind = FALSE)

# generate a PCA based of the DESeq object
plotPCA(vst, intgroup = c("treatment")) +
  aes(shape = treatment, color = treatment)

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


# Let's save these res_df
#write_rds(res_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230deseq_results_df.rds")

# shrunken log fold change results
#write_rds(res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230deseq_results_shrunken_lfc_df.rds")

res_df
res_shrunken_df

# pulling out control v 10 mM Zn results df (and all other comparisons)
cont_v_10_res_df <- res_df[res_df$result_name == "treatment_10_vs_0",]
cont_v_2.5_res_df <- res_df[res_df$result_name == "treatment_2.5_vs_0",]
cont_v_1_res_df <- res_df[res_df$result_name == "treatment_1_vs_0",]
cont_v_0.1_res_df <- res_df[res_df$result_name == "treatment_0.1_vs_0",]

# same for shrunken data
# pulling out control v 10 mM Zn results df (and all other comparisons)
cont_v_10_res_shrunken_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_10_vs_0",]
cont_v_2.5_res_shrunken_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_2.5_vs_0",]
cont_v_1_res_shrunken_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_1_vs_0",]
cont_v_0.1_res_shrunken_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_0.1_vs_0",]

# save each of these as a csv
#write.csv(cont_v_10_res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_cont_v_10_sh.csv")
#write.csv(cont_v_5_res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_cont_v_5_sh.csv")
#write.csv(cont_v_2.5_res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_cont_v_2_5_sh.csv")
#write.csv(cont_v_1_res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_cont_v_1_sh.csv")
#write.csv(cont_v_0.1_res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_cont_v_0_1_sh.csv")

# Let's make list of all SIG gnees with  P <0.05 and FC > 1.5 in the comparisons
cont_v_10_sig_genes <- cont_v_10_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_2.5_sig_genes <- cont_v_2.5_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_1_sig_genes <- cont_v_1_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_0.1_sig_genes <- cont_v_0.1_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)

sum(cont_v_10_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_10_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_2.5_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_2.5_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_1_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_1_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_0.1_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_0.1_sig_genes$log2FoldChange <= -1.5)

# shrunken
cont_v_10_shrunken_sig_genes <- cont_v_10_res_shrunken_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_2.5_shrunken_sig_genes <- cont_v_2.5_res_shrunken_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_1_shrunken_sig_genes <- cont_v_1_res_shrunken_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)
cont_v_0.1_shrunken_sig_genes <- cont_v_0.1_res_shrunken_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >1.5)

sum(cont_v_10_shrunken_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_10_shrunken_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_2.5_shrunken_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_2.5_shrunken_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_1_shrunken_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_1_shrunken_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_0.1_shrunken_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_0.1_shrunken_sig_genes$log2FoldChange <= -1.5)

# Now we filter res_shrunken for gene_id column and treatment
cont_v_10_genes_to_plot <- res_shrunken_df %>%
  filter(gene_id %in% cont_v_10_shrunken_sig_genes$gene_id, 
         result_name %in% c("treatment_0.1_vs_0", "treatment_1_vs_0", 
                            "treatment_2.5_vs_0", "treatment_10_vs_0"))

# We need a matrix for heatmap and converting genes/values from lines above to matrix
lfc_matrix <- cont_v_10_genes_to_plot %>% 
  dplyr::select(gene_id, log2FoldChange, result_name) %>% 
  pivot_wider(names_from = "result_name", values_from = "log2FoldChange") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

colnames(lfc_matrix)
# set order of columns to plot
lfc_matrix <- lfc_matrix[, c(1, 2, 4, 5, 3)]

# plot heat map of genes that are significantly up/dn regulated in 0 vs 10 - values for all treatments
pheatmap(lfc_matrix, show_rownames = F, breaks = seq(-3, 3, length.out = 100),
                   cluster_cols = F, border_color = NA)

# make a volcano plot of the differentiall expressed genes
keyvals <- ifelse(
  cont_v_10_res_shrunken_df$log2FoldChange < -1.5 & cont_v_10_res_shrunken_df$padj <= 0.05, 'blue',
  ifelse(cont_v_10_res_shrunken_df$log2FoldChange > 1.5 & cont_v_10_res_shrunken_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_10_res_shrunken_df,
                lab = rownames(cont_v_10_res_shrunken_df),
                selectLab = c('1'),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-10, 10), ylim=c(0,30))

# 2.5mM Zn
keyvals <- ifelse(
  cont_v_2.5_res_shrunken_df$log2FoldChange < -1.5 & cont_v_2.5_res_shrunken_df$padj <= 0.05, 'blue',
  ifelse(cont_v_2.5_res_shrunken_df$log2FoldChange > 1.5 & cont_v_2.5_res_shrunken_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_2.5_res_shrunken_df,
                lab = rownames(cont_v_2.5_res_shrunken_df),
                selectLab = c('1'),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-10, 10), ylim=c(0,30))

# 1mM Zn
keyvals <- ifelse(
  cont_v_1_res_shrunken_df$log2FoldChange < -1.5 & cont_v_1_res_shrunken_df$padj <= 0.05, 'blue',
  ifelse(cont_v_1_res_shrunken_df$log2FoldChange > 1.5 & cont_v_1_res_shrunken_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_1_res_shrunken_df,
                lab = rownames(cont_v_1_res_shrunken_df),
                selectLab = c(''),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-10, 10), ylim=c(0,30))


# 0.1mM Zn
keyvals <- ifelse(
  cont_v_0.1_res_shrunken_df$log2FoldChange < -1.5 & cont_v_0.1_res_shrunken_df$padj <= 0.05, 'blue',
  ifelse(cont_v_0.1_res_shrunken_df$log2FoldChange > 1.5 & cont_v_0.1_res_shrunken_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_0.1_res_shrunken_df,
                lab = rownames(cont_v_0.1_res_shrunken_df),
                selectLab = c(''),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-10, 10), ylim=c(0,30))

## GO enrichment analysis
## R script -> Stomentosus230_TopGO
