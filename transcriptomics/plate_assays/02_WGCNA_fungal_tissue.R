# WGCNA for S. tomentosus 230 isolate tissue ONLY in increasing Zn conditions

library(WGCNA)
library(gridExtra)
library(dplyr)
library(DESeq2)
library(ggplot2)

# counts sheet

ref230 <- read.table("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Transcriptome Analysis/counts.txt",
                     header = T, sep = '\t', fill = TRUE)

ref230_ <- ref230[,-1]
rownames(ref230_) <- ref230[,1]

# make a df for tissue 
countsall <- ref230_ %>%
  dplyr::select(grep("_tissue", colnames(ref230_), value = T))
colnames(countsall)

# remove outliers previously identified
exclude_samples <- c(
  "X220_tissue_0Zn_Rep1",
  "X230_tissue_0Zn_Rep4",
  "X220_tissue_1Zn_Rep1",
  "X230_tissue_1Zn_Rep3",
  "X230_tissue_5Zn_Rep1",
  "X230_tissue_5Zn_Rep2",
  "X230_tissue_5Zn_Rep3",
  "X230_tissue_5Zn_Rep4")

countsall <- countsall %>%
  dplyr::select(-all_of(exclude_samples))

colnames(countsall)
# hierarchical clustering 
htree <- hclust(dist(t(countsall)), method = "average")
plot(htree)
# indactes that there are no major outliers

# make col data
sample_names <- colnames(countsall)

# Function to extract isolate from the name
extract_info <- function(sample) {
  # Extract treatment (digits before 'Zn')
  treatment <- sub(".*_([0-9]+\\.?[0-9]*)Zn_.*", "\\1", sample)
  isolate <- ifelse(grepl("control_pine", sample), "none", sub(".*X(\\d+).*", "\\1", sample))
  
  return(c(sample, treatment, isolate))
}

# Apply the function to all sample names
col_data <- as.data.frame(t(sapply(sample_names, extract_info)))
colnames(col_data) <- c("", "treatment", "isolate")
row.names(col_data) <- NULL
col_data$treatment <- as.factor(col_data$treatment)
print(col_data)

#remove any NA datasets (change to 0s)
countsall[is.na(countsall)] <- 0
# create dds object (deseq)
dds <- DESeqDataSetFromMatrix(countData = countsall, 
                              colData = col_data,
                              design = ~1) # not specifying model

# remove all genes with <15 count in more than 75% of samples
n <- ncol(countsall)
cutoff <- n*0.75

# recommended by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) >= cutoff,]
number_genes <- nrow(dds75) # how many genes are left
number_genes

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) 

# needs to be transposed for WGCNA analysis
norm.counts <- norm.counts %>%
  t()

# WGCNA

# Network construction
# Choose a set of soft-thresholding powers
power <- c(c(1:12), seq(from = 14, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
# here y intercept is 0.8 and we want to choose a power that has an R2 of above 0.8 
# ie above the line (but not excessively large)

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
# we want a power with minimal mean connectivity

# look at both graphs at once to decide which power to choose
grid.arrange(a1, a2, nrow = 2)
# if no suitable = the authors recommend using a power of 18 
# for signed networks for a sample size between 
# 20 and 30 in case the scale free topology fit index fails to reach values
# from this we will choose power of 18 for our soft threshold

# but here - power 11 looks ideal

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 11 # from above
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = number_genes, # 10309 genes in our dds
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25, 
                          numericLabels = FALSE,
                          randomSeed = 1916,
                          verbose = 3)


cor <- temp_cor

moduleLabels <- bwnet$colors
moduleColors <- labels2colors(bwnet$colors)
MEs <- bwnet$MEs
geneTree <- bwnet$dendrograms[[1]]


#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, moduleLabels,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#gives table with numbers of genes in each module
# and plots it 
table(moduleLabels)
moduleColData <- as.data.frame(table(moduleLabels))

ggplot(data = moduleColData,
       aes(x = moduleLabels, y = Freq, 
           group = moduleLabels, fill = moduleLabels)) +
  geom_bar(stat = "identity", color = "black", show.legend = F) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("black", "blue", "brown",         
                               "green",
                               "greenyellow", "grey", 
                               "magenta", 
                               "pink", "purple", "red",
                               "salmon", "tan", 
                               "turquoise", "yellow"))

# Define numbers of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

col_data$treatment <- as.numeric(as.character(col_data$treatment))
col_data$isolate <- as.factor(col_data$isolate)
# make sensitive and tolerant a binary
# where tolerant = 1 and sentive = 0
col_data$tolerance <- ifelse(col_data$isolate == "230", 1, 0)

# select out only trait data
sampleDataZnforheatmap <- col_data %>% dplyr::select(treatment, tolerance)

## Recalculate MEs with color labels
MEs0 <- moduleEigengenes(norm.counts, moduleLabels)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, sampleDataZnforheatmap, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# will be changing margins for plots here so save original margins to return to them
old_par <- par(no.readonly = TRUE)
# The function “plotEigengeneNetworks” creates a heatmap of adjacencies 
# among eigengenes and a dendrogram for their relationship
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

# return margins to default
par(old_par)
## Will display correlations and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sampleDataZnforheatmap),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# modules of interest 
# sig negative correlation = pink? borderline
# sig positive correlation = purple

# Define variable Zn treatment containing the zn treatment column of the sample data file
Zn <- as.data.frame(sampleDataZnforheatmap$treatment)
names(Zn) <- "Zn"
geneTraitSignificanceZn <- as.data.frame(cor(norm.counts, Zn, use = "p"))
GSPvalueZn <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceZn), nSamples))
names(geneTraitSignificanceZn) <- paste("GS.", names(Zn), sep="")
names(GSPvalueZn) <- paste("p.GS.", names(Zn), sep="")
head(GSPvalueZn)

# do same for tolerance
Tol <- as.data.frame(sampleDataZnforheatmap$tolerance)
names(Tol) <- "Tol"
geneTraitSignificanceTol <- as.data.frame(cor(norm.counts, Tol, use = "p"))
GSPvalueTol <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceTol), nSamples))
names(geneTraitSignificanceTol) <- paste("GS.", names(Tol), sep="")
names(GSPvalueTol) <- paste("p.GS.", names(Tol), sep="")
head(GSPvalueTol)

#Calculate the module membership and the associated p-values
geneModuleMembership <- as.data.frame(cor(norm.counts, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
modNames <- substring(names(MEs), 3) #extract module names
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")


# Make a graph of module membership vs trait p value to highlight important genes

# modules correlated with Zn treatment
# sig pos corr - purple, greenyellow, green
# sig neg corr - pink

module <- "purple" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceZn[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Zn treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "greenyellow"
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceZn[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Zn treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "green" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceZn[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Zn treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "pink" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceZn[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Zn treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# modules correlated with tolerance - change the geneTraitSignificanceZn to Tol
# sig pos corr - blue, greenyellow, magenta
# sig neg corr - yellow, grey (pink is borderline)

module <- "blue" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "greenyellow" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "magenta" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "yellow" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "grey" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "pink" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module <- "brown" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceTol[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Tolerance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# output lists of genes in each module
genes.pvals.Zn <- tibble::rownames_to_column(GSPvalueZn, "gene")
genes.MM.Zn <- tibble::rownames_to_column(GSPvalueZn, "gene")
# transform p values -log10 = higher = more significant
genes.pvals.Zn$rank <- -log10(genes.pvals.Zn$p.GS.Zn)

genes.pvals.Tol <- tibble::rownames_to_column(GSPvalueTol, "gene")
genes.MM.ZTol <- tibble::rownames_to_column(GSPvalueTol, "gene")
# transform p values -log10 = higher = more significant
genes.pvals.Tol$rank <- -log10(genes.pvals.Tol$p.GS.Tol)


# output the lists and change the gene col to be just protein ID
# sig pos corr Zn treatment - purple, greenyellow, green
# sig neg corr Zn treatment - pink

purpleZnrank <- genes.pvals.Zn[moduleLabels=="purple",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

greenyellowZnrank <- genes.pvals.Zn[moduleLabels=="greenyellow",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

greenZnrank <- genes.pvals.Zn[moduleLabels=="green",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

pinkZnrank <- genes.pvals.Zn[moduleLabels=="pink",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

# same but for tolerance modules
# sig pos corr - blue, greenyellow, magenta
# sig neg corr - yellow, grey (pink is borderline)
blueTolrank <- genes.pvals.Tol[moduleLabels=="blue",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

greenyellowTolrank <- genes.pvals.Tol[moduleLabels=="greenyellow",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

magentaTolrank <- genes.pvals.Tol[moduleLabels=="magenta",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

yellowTolrank <- genes.pvals.Tol[moduleLabels=="yellow",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

greyTolrank <- genes.pvals.Tol[moduleLabels=="grey",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

pinkTolrank <- genes.pvals.Tol[moduleLabels=="pink",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))

brownTolrank <- genes.pvals.Tol[moduleLabels=="brown",] %>%
  arrange( ,rank, desc(rank)) %>%
  mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))


# output the rank lists as csv - if wanted
# write.csv(purpleZnrank, "~/Desktop/WGCNA_tissue/230/purpleZnrank.csv") 


### GO ANALYSIS ON MODULES
library(topGO)
# read the GO mapping file - contains which GO terms relate to which genes
gene_to_go_mapping <- readMappings(file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Stomentosus_gene_to_go.txt", 
                                   sep = "\t", IDsep = ",")

# all genes are all the genes with GO annotations
allGenes <- names(gene_to_go_mapping)

# geneList is a vector with the gene IDs as names
# and a 1/0 if genes are present in the rank list

# modules correlated with Zn treatment
# sig pos corr - purple, greenyellow, green
# sig neg corr - pink
# modules correlated with tolerance - ** change to -Tolrank
# sig pos corr - blue, greenyellow, magenta
# sig neg corr - yellow, grey (pink is borderline)

# doesnt matter the rank here as only looking for presence absence in the module


moduleGO <- brownTolrank
geneList <- factor(as.integer(allGenes %in% moduleGO$gene))
names(geneList) <- allGenes

# create topGO object - Molecular Function
MFdata <- new("topGOdata",
              ontology = "MF",  
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene_to_go_mapping)

resultTableMF <- runTest(MFdata, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(MFdata, 
                     weightFisher = resultTableMF, 
                     topNodes = 20)
allResMF

# Biological process
BPdata <- new("topGOdata",
              ontology = "BP",  
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene_to_go_mapping)

resultTableBP <- runTest(BPdata, algorithm = "weight01", statistic = "fisher")
allResBP <- GenTable(BPdata, 
                     weightFisher = resultTableBP, 
                     topNodes = 20)
allResBP

# Cellular component
CCdata <- new("topGOdata",
              ontology = "CC",  
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene_to_go_mapping)

resultTableCC <- runTest(CCdata, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(CCdata, 
                     weightFisher = resultTableCC, 
                     topNodes = 20)
allResCC

##
##
##
# Intramodular analysis

# GSEA
# have list of genes significantly correlated with Zn treatment from above
# use this to create a pre-ranked list based on -log10 p val to to GSEA
forGSEAZn <- tibble::rownames_to_column(GSPvalueZn, "gene")
colnames(forGSEAZn) <- c("gene","p.val")
forGSEAZn <- forGSEAZn %>% mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))
forGSEAZn$rank <- -log10(forGSEAZn$p.val)
forGSEAZn <- arrange(forGSEAZn, desc(rank))
write.csv(forGSEAZn, "~/Desktop/Stomentosus_Zn_Transcriptomics/WGCNA_tissue/Zn_treatment_correlated_genes.csv") 

# GO analysis of intramodulara analysis
# first, get only the name and the p value
Genes_Universe <- forGSEAZn %>%
  dplyr::select(gene, p.val)

# sort by the membership value
Genes_Universe <- Genes_Universe[order(Genes_Universe$p.val), ]

# convert to topGO's genelist format
topgoUniverse <- Genes_Universe$p.val
names(topgoUniverse) <- Genes_Universe$gene
topgoUniverse <- na.omit(topgoUniverse)

# Gene Selection needs to run off a function 
# This function returns the genes with p value under 0.05
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)}

my_go_dataBP <- new("topGOdata",
                    ontology    = "BP",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableBP <- runTest(my_go_dataBP, algorithm = "weight01", statistic = "fisher")
allResBP <- GenTable(my_go_dataBP, 
                     weightFisher = resultTableBP, 
                     topNodes = 20)
allResBP

my_go_dataMF <- new("topGOdata",
                    ontology    = "MF",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableMF <- runTest(my_go_dataMF, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(my_go_dataMF, 
                     weightFisher = resultTableMF, 
                     topNodes = 20)
allResMF

my_go_dataCC <- new("topGOdata",
                    ontology    = "CC",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableCC <- runTest(my_go_dataCC, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(my_go_dataCC, 
                     weightFisher = resultTableCC, 
                     topNodes = 20)
allResCC

##
##
## Now for tolerance correlated genes
GSPvalueTol %>% 
  head(25)

forGSEATol <- as.data.frame(GSPvalueTol)
forGSEATol <- tibble::rownames_to_column(forGSEATol, "gene")
colnames(forGSEATol) <- c("gene","p.val")
forGSEATol <- forGSEATol %>% mutate(gene = sub(".*proteinId=([0-9]+).*", "\\1", gene))
forGSEATol$rank <- -log10(forGSEATol$p.val)
forGSEATol <- arrange(forGSEATol, desc(rank))
write.csv(forGSEATol, "~/Desktop/Stomentosus_Zn_Transcriptomics/WGCNA_tissue/tolerance_correlated_genes.csv") 

# perform GSEA preranked manually on GSEA GUI 

Genes_Universe <- forGSEATol %>%
  dplyr::select(gene, p.val)

# sort by the membership value
Genes_Universe <- Genes_Universe[order(Genes_Universe$p.val), ]

# convert to topGO's genelist format
topgoUniverse <- Genes_Universe$p.val
names(topgoUniverse) <- Genes_Universe$gene
topgoUniverse <- na.omit(topgoUniverse)

# Gene Selection needs to run off a function 
# This function returns the genes with p value under 0.05
topDiffGenes <- function(allScore) {
  return(allScore < 0.05)}

my_go_dataBP <- new("topGOdata",
                    ontology    = "BP",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableBP <- runTest(my_go_dataBP, algorithm = "weight01", statistic = "fisher")
allResBP <- GenTable(my_go_dataBP, 
                     weightFisher = resultTableBP, 
                     topNodes = 20)
allResBP

my_go_dataMF <- new("topGOdata",
                    ontology    = "MF",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableMF <- runTest(my_go_dataMF, algorithm = "weight01", statistic = "fisher")
allResMF <- GenTable(my_go_dataMF, 
                     weightFisher = resultTableMF, 
                     topNodes = 20)
allResMF

my_go_dataCC <- new("topGOdata",
                    ontology    = "CC",
                    geneSel     = topDiffGenes, # this needs to be a function
                    allGenes    = topgoUniverse, # this is all the genes in the expt
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 10) # Modify to reduce/increase stringency.

resultTableCC <- runTest(my_go_dataCC, algorithm = "weight01", statistic = "fisher")
allResCC <- GenTable(my_go_dataCC, 
                     weightFisher = resultTableCC, 
                     topNodes = 20)
allResCC
