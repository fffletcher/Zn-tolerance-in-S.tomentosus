# Stomentosus DE gene Expression
# using reads aligned to 230 genome (better genome)

# load libraries
library(tidyverse)
library(dplyr)
library(pheatmap)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)

# comparison between 220 (sensitive) and 230 (tolerant) isoaltes
# in 0 mM Zn conditions - constituitive tolerance? 

# First look at tissue - look at 230
# counts sheet

ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# make a df for tissue from 230
counts <- ref230_ %>%
  dplyr::select(grep("tissue_0Zn", colnames(ref230_), value = T))
colnames(counts)


# remove outliers - identified previously
# 220 outliers = 0Zn_Rep1
# 230 outliers = 0Zn_Rep4
tissuecounts <- counts %>%
  dplyr::select(!grep("X220_tissue_0Zn_Rep1", colnames(counts), value = T)) %>%
  dplyr::select(!grep("X230_tissue_0Zn_Rep4", colnames(counts), value = T))
colnames(tissuecounts)

# make col data files from samples in the df
sample_names <- colnames(tissuecounts)

# Function to extract isolate from the name
extract_info <- function(sample) {
  # Extract isolate (digits after 'X', or 'none' if 'control_pine' is in the name)
  isolate <- ifelse(grepl("control_pine", sample), "none", sub(".*X(\\d+).*", "\\1", sample))
  
  return(c(sample, isolate))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("","isolate")
row.names(col_data) <- NULL
col_data$isolate <- as.factor(col_data$isolate)
print(col_data)

#remove any NA datasets (change to 0s)
tissuecounts[is.na(tissuecounts)] <- 0
# Create DESeq object (note design only by treatment as they are all the same isolate)
ddsMat <- DESeqDataSetFromMatrix(countData = tissuecounts,
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
  aes(shape = isolate, color = isolate)

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

res_df
res_shrunken_df

# Let's save these res_df
#write_rds(res_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/230vs230_0Zn_deseq_results_df.rds")

# shrunken log fold change results
# write_rds(res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/230vs230_0Zn_results_shrunken_lfc_df.rds")

# save the shrunken as a csv too
write.csv(res_shrunken_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/230vs230_0Zn_sh.csv")

# how many genes up and down regualted
sum(res_shrunken_df$log2FoldChange >= 1.5)
sum(res_shrunken_df$log2FoldChange <= -1.5)

# volcano plot
keyvals <- ifelse(
  res_shrunken_df$log2FoldChange < -1.5 & res_shrunken_df$padj <= 0.05, '#E69F00',
  ifelse(res_shrunken_df$log2FoldChange > 1.5 & res_shrunken_df$padj <= 0.05, '#0072B2',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#0072B2'] <- 'upregulated'
names(keyvals)[keyvals == '#E69F00'] <- 'down-regualted'

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
                gridlines.minor = F) +
  ggplot2::coord_cartesian(xlim=c(-20, 20), ylim=c(0,40)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  annotate("text", x = -13, y = 39, label = "Genes more highly expressed in ZnS", 
           color = "#E69F00", size = 3, fontface = "bold") +
  annotate("text", x = 13, y = 39, label = "Genes more highly expressed in ZnT", 
           color = "#0072B2", size = 3, fontface = "bold")



# output a list of genes ordered by log2fc - lowest first
# for use in GSEA
# no headers and only protin ID
# save as .txt file with the .rnk file type
# extract protein ID
df <- res_shrunken_df %>%
  mutate(gene_id = sub(".*proteinId=([0-9]+).*", "\\1", gene_id))

rank <- df %>%
  dplyr::select(gene_id, log2FoldChange) %>%
  arrange(log2FoldChange)
rank$gene_id <- as.numeric(rank$gene_id)
names(rank) <- NULL

write.csv(rank, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/220v230.csv")

