# Pine DE gene Expression - BIOASSAY
# fungal inoculation versus no inocuation


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

ref <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref_ <- ref[,-1]
rownames(ref_) <- ref[,1]

# dont include contol pine samples - we want to look only at Zn treatment effect on inoculated pines
# there was no real difference in isolate so we can look at them together

pine <- ref_ %>%
  dplyr::select(grep("_0Zn", colnames(ref_), value = T))
pine <- pine %>%
  dplyr::select(!grep("tissue", colnames(pine), value = T))

# exclude 230 rep 4 and 220 rep 1 - outliers
pine <- pine %>%
  dplyr::select(!grep("X230_bioassay_0Zn_Rep4", colnames(pine), value = T)) 
pine <- pine %>%
  dplyr::select(!grep("X220_bioassay_0Zn_Rep1", colnames(pine), value = T))

colnames(pine)

# make col data files from samples in the df
sample_names <- colnames(pine)

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
pine[is.na(pine)] <- 0

# Create DESeq object (note design only by treatment)
ddsMat <- DESeqDataSetFromMatrix(countData = pine,
                                 colData = col_data,
                                 design = ~ inocualted)

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
  #geom_text(label = colnames(vst))

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
  res_shrunken_df$log2FoldChange <= -1.5 & res_shrunken_df$padj < 0.05, '#1b9f77',
  ifelse(res_shrunken_df$log2FoldChange >= 1.5 & res_shrunken_df$padj < 0.05, '#D55E00',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#D55E00'] <- 'upregulated'
names(keyvals)[keyvals == '#1b9f77'] <- 'down-reguated'

EnhancedVolcano(res_shrunken_df,
                lab = rownames(res_shrunken_df),
                selectLab = c('1'),
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
gene_to_go_mapping <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Pinuscontorta_gene_to_go.txt", 
                                   sep = "\t", IDsep = ",")


names(gene_to_go_mapping) <- gsub('"', '', names(gene_to_go_mapping))

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
sig_go_df
#write.csv(sig_go_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_DEG_GO_10mM.csv")


# what about significant differences between inoculation with ZnT and ZnS
fungus <- ref_ %>%
  dplyr::select(grep("_bioassay_0Zn", colnames(ref_), value = T))

# exclude 230 rep 4 and 220 rep 1
fungus <- fungus %>%
  dplyr::select(!grep("X230_bioassay_0Zn_Rep4", colnames(fungus), value = T)) 
fungus <- fungus %>%
  dplyr::select(!grep("X220_bioassay_0Zn_Rep1", colnames(fungus), value = T))

colnames(fungus)

# make col data files from samples in the df
sample_names <- colnames(fungus)

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
fungus[is.na(fungus)] <- 0

# Create DESeq object (note design only by treatment)
ddsMat <- DESeqDataSetFromMatrix(countData = fungus,
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
#geom_text(label = colnames(vst))

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
  ggplot2::coord_cartesian(xlim=c(-15, 15))

# GO enrichment
library(topGO)
# read the GO mapping file - with IDsep specification
gene_to_go_mapping <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Pinuscontorta_gene_to_go.txt", 
                                   sep = "\t", IDsep = ",")


names(gene_to_go_mapping) <- gsub('"', '', names(gene_to_go_mapping))

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
sig_go_df
#write.csv(sig_go_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/230_DEG_GO_10mM.csv")