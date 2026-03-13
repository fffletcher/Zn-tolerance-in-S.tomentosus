## Making counts files for S. tomentosus transcriptomics
# Want each isolate in each set seperate - 220 and 230
library(DESeq2)
library(dplyr)
library(ggplot2)

## load in fungal reads aligned to 230 reference - includes bioassay and tissue
ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# make a df for bioassay only 
tissue_ref230 <- ref230_ %>%
  dplyr::select(grep("bioassay", colnames(ref230_), value = T), 
                grep("control", colnames(ref230_), value = T))

# make col data file
sample_names <- colnames(tissue_ref230)

# Function to extract isolate and treatment from the name
extract_info <- function(sample) {
  # Extract isolate (digits after 'X', or 'none' if 'control_pine' is in the name)
  isolate <- ifelse(grepl("control_pine", sample), "none", sub(".*X(\\d+).*", "\\1", sample))
  
  # Extract treatment (digits before 'Zn')
  treatment <- sub(".*_([0-9]+\\.?[0-9]*)Zn_.*", "\\1", sample)
  
  return(c(sample, isolate, treatment))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("","isolate", "treatment")
row.names(col_data) <- NULL
col_data$isolate <- as.factor(col_data$isolate)
col_data$treatment <- factor(col_data$treatment,
                             levels = c("0", "0.1", "1", "2.5", "5", "10"))

print(col_data)


# PCA

#remove any NA datasets (change to 0s)
tissue_ref230[is.na(tissue_ref230)] <- 0
# Create DESeq object (note design only by treatment as they are all the same isolate)
ddsMat <- DESeqDataSetFromMatrix(countData = tissue_ref230,
                                 colData = col_data,
                                 design = ~ treatment + isolate)

# Filtering: removing rows of the DESeqDataSet that have no counts, 
# or only a single count across all samples.
keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

# rlog transformation
# variance stabilizing transformation (VST)
vst <- vst(ddsMat, blind = FALSE)
plotPCA(vst, intgroup = c("treatment", "isolate")) +
  aes(shape = treatment, color = isolate)

plotPCA(vst, intgroup = c("treatment", "isolate")) +
  aes(shape = treatment, color = isolate) +
  scale_color_manual(
    values = c("none" = "grey" ,"220" = "#E69F00", "230" = "#0072B2"),
    labels = c("220" = "ZnS", "230" = "ZnT"),  # your desired mapping
    name = NULL
  ) +
  scale_shape_manual(
    values = c(8, 15, 16, 18),
    name = "Zn treatment (mM)"
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  geom_point(size = 4) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#E69F00", "#0072B2", "darkgrey"),
                    labels = c("none" = "None",
                               "220" = "ZnS",
                               "230" = "ZnT"),
                    name = "Fungus")

# PERMANOVA 
library(factoextra)
library(vegan)

test <- t(assay(vst))
test <- as.data.frame(test)
cols <- ncol(test)

# add sample and treatment columns
test$sample <- rownames(test)
col_data1 <- col_data
names(col_data1)[1] <- "sample_id" 
test$treatment <- col_data1$treatment[match(test$sample, col_data1$sample_id)]
test$isolate <- col_data1$isolate[match(test$sample, col_data1$sample_id)]

#princopal components analysis
scaled_test <- prcomp(test[c(1:cols)], scale = TRUE, center = TRUE)
# outputs eigenvalues
get_eig(scaled_test)

# scale data
vegan <- scale(test[c(1:cols)])

# permanova
# use "eu" method as we scaled our data
adonis2(vegan ~ isolate, data = test, method = 'eu') # only isolate
adonis2(vegan ~ treatment, data = test, method = 'eu') # only treatment

# interaction
adonis2(vegan ~ isolate * treatment, data = test, method = 'eu') 


# pairwise
library(pairwiseAdonis)
pairwise.adonis2(vegan ~ treatment, data = test, method = 'eu')
pairwise.adonis2(vegan ~ isolate, data = test, method = 'eu')

# plot heatmap for all sample to sample distances
library(pheatmap)
library(RColorBrewer)

gsampleDists <- dist(t(assay(vst)))
gsampleDistMatrix <- as.matrix(gsampleDists)
rownames(gsampleDistMatrix) <- colnames(vst)
colnames(gsampleDistMatrix) <- NULL

colors <- colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)

pheatmap(gsampleDistMatrix,
         clustering_distance_rows = gsampleDists,
         clustering_distance_cols = gsampleDists,
         color = colors)

# seperate PCA for each isolate
# only 230 isolate
i230_ref230 <- tissue_ref230 %>%
  dplyr::select(grep("X230", colnames(tissue_ref230), value = T))
col_data230 <- col_data[col_data$isolate=="230",]
# only 220 isolate
i220_ref230 <- tissue_ref230 %>%
  dplyr::select(grep("X220", colnames(tissue_ref230), value = T))
col_data220 <- col_data[col_data$isolate=="220",]


ddsMat230 <- DESeqDataSetFromMatrix(countData = i230_ref230,
                                    colData = col_data230,
                                    design = ~ treatment)

ddsMat220 <- DESeqDataSetFromMatrix(countData = i220_ref230,
                                    colData = col_data220,
                                    design = ~ treatment)

keep230 <- rowSums(counts(ddsMat230)) > 1
ddsMat230 <- ddsMat230[keep230,]

keep220 <- rowSums(counts(ddsMat220)) > 1
ddsMat220 <- ddsMat220[keep220,]

vst230 <- vst(ddsMat230, blind = FALSE)
plotPCA(vst230, intgroup = c("treatment")) +
  aes()

vst220 <- vst(ddsMat220, blind = FALSE)
plotPCA(vst220, intgroup = c("treatment")) +
  aes()

# PERMANOVA 
test220 <- t(assay(vst220))
test220 <- as.data.frame(test220)
cols220 <- ncol(test220)

test230 <- t(assay(vst230))
test230 <- as.data.frame(test230)
cols230 <- ncol(test230)

# add sample and treatment columns
test220$sample <- rownames(test220)
col_data220_1 <- col_data220
names(col_data220_1)[1] <- "sample_id" 
test220$treatment <- col_data220_1$treatment[match(test220$sample, col_data220_1$sample_id)]
test220$isolate <- col_data220_1$isolate[match(test220$sample, col_data220_1$sample_id)]

test230$sample <- rownames(test230)
col_data230_1 <- col_data230
names(col_data230_1)[1] <- "sample_id" 
test230$treatment <- col_data230_1$treatment[match(test230$sample, col_data230_1$sample_id)]
test230$isolate <- col_data230_1$isolate[match(test230$sample, col_data230_1$sample_id)]

# scale and get eigenvalues
scaled_test220 <- prcomp(test220[c(1:cols220)], scale = TRUE, center = TRUE)
get_eig(scaled_test220)
vegan220 <- scale(test220[c(1:cols220)])

scaled_test230 <- prcomp(test230[c(1:cols230)], scale = TRUE, center = TRUE)
get_eig(scaled_test230)
vegan230 <- scale(test230[c(1:cols230)])

# permanova
# use "eu" method as we scaled our data
adonis2(vegan220 ~ treatment, data = test220, method = 'eu') # only treatment (it is the only variable)
adonis2(vegan230 ~ treatment, data = test230, method = 'eu')

# pairwise
library(pairwiseAdonis)
pairwise.adonis2(vegan220 ~ treatment, data = test220, method = 'eu')
pairwise.adonis2(vegan230 ~ treatment, data = test230, method = 'eu')

# plot heatmap for all sample to sample distances
gsampleDists220 <- dist(t(assay(vst220)))
gsampleDistMatrix220 <- as.matrix(gsampleDists220)
rownames(gsampleDistMatrix220) <- colnames(vst220)
colnames(gsampleDistMatrix220) <- NULL

gsampleDists230 <- dist(t(assay(vst230)))
gsampleDistMatrix230 <- as.matrix(gsampleDists230)
rownames(gsampleDistMatrix230) <- colnames(vst230)
colnames(gsampleDistMatrix230) <- NULL

pheatmap(gsampleDistMatrix220,
         clustering_distance_rows = gsampleDists220,
         clustering_distance_cols = gsampleDists220,
         color = colors)

pheatmap(gsampleDistMatrix230,
         clustering_distance_rows = gsampleDists230,
         clustering_distance_cols = gsampleDists230,
         color = colors)

