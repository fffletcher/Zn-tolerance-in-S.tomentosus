# Stomentosus DE gene Expression - BIOASSAY
# using reads aligned to FUNGAL 230 genome (better genome)
# What comparisons are wanted? 

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

# First look at tissue - look at 230
# counts sheet

ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# dont include contol pine samples - we want to look only at Zn treatment effect on inoculated pines
# there was no real difference in isolate so we can look at them together

tissue_ref230 <- ref230_ %>%
  dplyr::select(grep("bioassay", colnames(ref230_), value = T))

colnames(tissue_ref230)

# make col data files from samples in the df
sample_names <- colnames(tissue_ref230)

# Function to extract isolate from the name
extract_info <- function(sample) {
  # Extract treatment (digits before 'Zn')
  treatment <- sub(".*_([0-9]+\\.?[0-9]*)Zn_.*", "\\1", sample)
  isolate <- ifelse(grepl("control_pine", sample), "none", sub(".*X(\\d+).*", "\\1", sample))
  
  return(c(sample, treatment, isolate))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("sample", "treatment", "isolate")
row.names(col_data) <- NULL
col_data$treatment <- as.factor(col_data$treatment)
print(col_data)

#remove any NA datasets (change to 0s)
tissue_ref230[is.na(tissue_ref230)] <- 0

# Create DESeq object (note design only by treatment)
ddsMat <- DESeqDataSetFromMatrix(countData = tissue_ref230,
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

# pulling out control v 10 mM Zn results df (and all other comparisons)
# use shrunken
cont_v_10_res_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_10_vs_0",]
cont_v_5_res_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_5_vs_0",]
cont_v_1_res_df <- res_shrunken_df[res_shrunken_df$result_name == "treatment_1_vs_0",]

# write these to csv and save for further anaysis
write.csv(cont_v_10_res_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/pine_fungal_inoculated_only_cont_v_10.csv")
write.csv(cont_v_5_res_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/pine_fungal_inoculated_only_cont_v_5.csv")
write.csv(cont_v_1_res_df, "~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/pine_fungal_inoculated_only_cont_v_1.csv")

# Let's make list of all SIG gnees with  P <0.05 and FC > 1.5 in the comparisons
cont_v_10_sig_genes <- cont_v_10_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.5)
cont_v_5_sig_genes <- cont_v_5_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.5)
cont_v_1_sig_genes <- cont_v_1_res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 1.5)

sum(cont_v_10_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_10_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_5_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_5_sig_genes$log2FoldChange <= -1.5)

sum(cont_v_1_sig_genes$log2FoldChange >= 1.5)
sum(cont_v_1_sig_genes$log2FoldChange <= -1.5)


# Now we filter res_shrunken for gene_id column and treatment
cont_v_10_genes_to_plot <- res_shrunken_df %>%
  filter(gene_id %in% cont_v_10_sig_genes$gene_id, 
         result_name %in% c("treatment_1_vs_0", 
                            "treatment_5_vs_0", 
                            "treatment_10_vs_0"))

# We need a matrix for heatmap and converting genes/values from lines above to matrix
lfc_matrix <- cont_v_10_genes_to_plot %>% 
  dplyr::select(gene_id, log2FoldChange, result_name) %>% 
  pivot_wider(names_from = "result_name", values_from = "log2FoldChange") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

colnames(lfc_matrix)
# set order of columns to plot
lfc_matrix <- lfc_matrix[, c(1, 3, 2)]

# plot heat map of genes that are significantly up/dn regulated in 0 vs 10 - values for all treatments
pheatmap(lfc_matrix, show_rownames = F, breaks = seq(-3, 3, length.out = 100),
         cluster_cols = F)

# make a volcano plot of the differentiall expressed genes
keyvals <- ifelse(
  cont_v_10_res_df$log2FoldChange <= -1.5 & cont_v_10_res_df$padj < 0.05, 'blue',
  ifelse(cont_v_10_res_df$log2FoldChange >= 1.5 & cont_v_10_res_df$padj < 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_10_res_df,
                lab = rownames(cont_v_10_res_df),
                selectLab = c('1'),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-15, 15))

# 5mM Zn
keyvals <- ifelse(
  cont_v_5_res_df$log2FoldChange <= -1.5 & cont_v_5_res_df$padj <= 0.05, 'blue',
  ifelse(cont_v_5_res_df$log2FoldChange >= 1.5 & cont_v_5_res_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_5_res_df,
                lab = rownames(cont_v_5_res_df),
                selectLab = c('1'),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-15, 15))

# 1mM Zn
keyvals <- ifelse(
  cont_v_1_res_df$log2FoldChange <= -1.5 & cont_v_1_res_df$padj <= 0.05, 'blue',
  ifelse(cont_v_1_res_df$log2FoldChange >= 1.5 & cont_v_1_res_df$padj <= 0.05, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'upregulated'
names(keyvals)[keyvals == 'blue'] <- 'down-regualted'

EnhancedVolcano(cont_v_1_res_df,
                lab = rownames(cont_v_1_res_df),
                selectLab = c('1'),
                x = 'log2FoldChange',
                y = 'padj',
                colCustom = keyvals,
                pCutoff = 0.05,
                FCcutoff = 1.5,
                colAlpha = 0.9,
                gridlines.major = F,
                gridlines.minor = T) +
  ggplot2::coord_cartesian(xlim=c(-15, 15))

# bar plots of degs across treatments
deg_data <- res_shrunken_df %>%
  mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1.5 ~ "up",
      padj < 0.05 & log2FoldChange < -1.5 ~ "down",
      TRUE ~ "not_significant"
    )
  )

# Set your desired x-axis order
deg_data$result_name <- factor(deg_data$result_name, 
                               levels = c("treatment_1_vs_0",
                                          "treatment_5_vs_0",
                                          "treatment_10_vs_0"))

# Count up/down regulated genes per condition
reg_counts <- deg_data %>%
  filter(regulation != "not_significant") %>%
  group_by(result_name, regulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(regulation == "down", -count, count))

# Scale counts to match y-axis range
ymin <- -15
ymax <- 15
max_count <- max(abs(reg_counts$count))
scale_factor <- ymax / max_count
reg_counts <- reg_counts %>%
  mutate(scaled_count = count * scale_factor)

# Separate significant and non-significant data
sig_data <- deg_data %>% filter(regulation != "not_significant")
nonsig_data <- deg_data %>% filter(regulation == "not_significant")

# Plot
p <- ggplot() +
  # Non-significant points
  geom_jitter(data = nonsig_data, aes(x = result_name, y = log2FoldChange, size = -log10(padj)),
              color = "grey", width = 0.2, alpha = 0.2) +
  # Significant points
  geom_jitter(data = sig_data, aes(x = result_name, y = log2FoldChange, color = regulation, size = -log10(padj)),
              width = 0.2, alpha = 0.8) +
  # DEG count bars
  geom_bar(data = reg_counts, aes(x = result_name, y = scaled_count, fill = regulation),
           stat = "identity", position = "identity", alpha = 0.3, width = 0.6) +
  # DEG count labels
  geom_text(
    data = reg_counts,
    aes(x = result_name, y = scaled_count, label = abs(count)),
    vjust = ifelse(reg_counts$count > 0, -1, 1),  # Above for up, below for down
    hjust = 0.5,  # Centered horizontally
    size = 3,
    inherit.aes = FALSE
  ) +
  # Scales and themes
  scale_color_manual(values = c("up" = "red", "down" = "blue")) +
  scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of Significant Genes")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Condition",
    y = "Log2 Fold Change",
    color = "Regulation",
    size = "-log10(padj)",
    fill = "Regulation"
  )

p

# make sure all x values are present even if no data
# and ensure same p val scale across isoaltes
p + scale_x_discrete(limits = c(
  #"treatment_0.1_vs_0", 
  "treatment_1_vs_0",
  #"treatment_2.5_vs_0",
  "treatment_5_vs_0",
  "treatment_10_vs_0")) +
  
  scale_size_continuous(
    name = "-log10(padj)",
    limits = c(0, 20),         # Set the same range across datasets
    breaks = c(5, 10, 15, 20), # Consistent legend ticks
    range = c(1, 5)            # Visual size range of points
  )

## venn diagram of up/down regualted genes - which are shared

# Define significance thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 1.5

# Get list of treatments
treatments <- unique(res_shrunken_df$result_name)

# Initialize lists for up/down regulated genes
upregulated_genes <- list()
downregulated_genes <- list()

# Loop through treatments and collect significant genes
for (treatment in treatments) {
  filtered <- res_shrunken_df %>%
    filter(result_name == treatment, padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  up_genes <- filtered %>%
    filter(log2FoldChange > lfc_cutoff) %>%
    pull(gene_id)
  
  down_genes <- filtered %>%
    filter(log2FoldChange < -lfc_cutoff) %>%
    pull(gene_id)
  
  if (length(up_genes) > 0) {
    upregulated_genes[[treatment]] <- up_genes
  }
  if (length(down_genes) > 0) {
    downregulated_genes[[treatment]] <- down_genes
  }
}


# Plot upregulated genes Venn diagram
if (length(upregulated_genes) >= 2) {
  p_up <- ggvenn(upregulated_genes,
                 fill_color = c("#E41A1C", "#FF7F00", "#4DAF4A", "#377EB8", "#984EA3"),
                 stroke_size = 0.5,
                 set_name_size = 5,
                 text_size = 4,
                 show_percentage = FALSE) +
    ggtitle("Upregulated Genes (padj < 0.05, LFC > 1.5)")
  print(p_up)
} else {
  message("Not enough treatments with upregulated genes to draw a Venn diagram.")
}

# Plot downregulated genes Venn diagram
if (length(downregulated_genes) >= 2) {
  p_down <- ggvenn(downregulated_genes,
                   fill_color = c("#377EB8", "#984EA3", "#A65628", "#999999", "#66C2A5"),
                   stroke_size = 0.5,
                   set_name_size = 5,
                   text_size = 4,
                   show_percentage = FALSE) +
    ggtitle("Downregulated Genes (padj < 0.05, LFC < -1.5)")
  print(p_down)
} else {
  message("Not enough treatments with downregulated genes to draw a Venn diagram.")
}


# Function to extract Venn regions from a list of sets
get_venn_regions <- function(gene_sets) {
  treatments <- names(gene_sets)
  combinations <- unlist(lapply(1:length(treatments), function(i) combn(treatments, i, simplify = FALSE)), recursive = FALSE)
  
  region_list <- list()
  
  for (combo in combinations) {
    # Genes in all treatments in combo
    shared_genes <- Reduce(intersect, gene_sets[combo])
    
    # Genes not in any other treatment
    other_treatments <- setdiff(treatments, combo)
    if (length(other_treatments) > 0) {
      other_genes <- Reduce(union, gene_sets[other_treatments])
      shared_genes <- setdiff(shared_genes, other_genes)
    }
    
    region_name <- paste(combo, collapse = " & ")
    region_list[[region_name]] <- shared_genes
  }
  
  return(region_list)
}

# Get Venn regions for upregulated and downregulated genes
up_regions <- get_venn_regions(upregulated_genes)
down_regions <- get_venn_regions(downregulated_genes)

# View results
up_regions
down_regions

## GO enrichment analysis
## R script -> Stomentosus230_fungalBioassay_TopGO
