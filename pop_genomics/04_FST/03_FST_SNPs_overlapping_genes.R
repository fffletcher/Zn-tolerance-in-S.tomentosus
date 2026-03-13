# Load libraries
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

setwd("~/Desktop/Stomentosus_PopGen/FST")

# check if SNPs in the windows overlap with genes or promotors

# first for FST = 1 windows only
# Read your FST regions
snps_df <- read.csv("snps_in_fst_one_windows.csv")  

# Read your GFF/GTF annotation
gff <- import("~/Desktop/Stomentosus_PopGen/Suillustomentosus.gff")  

# Convert SNPs to GRanges
snp_gr <- GRanges(seqnames = snps_df$CHROM,
                  ranges = IRanges(start = snps_df$POS, end = snps_df$POS),
                  SNP = snps_df$SNP)

# Filter gene features
gene_gr <- gff[gff$type == "gene"]

# Define promotor regions (1000 bp upstream of gene start, strand-aware)
promotor_gr <- promoters(gene_gr, upstream = 1000, downstream = 0)

# Initialize annotation column
snps_df$Annotation <- "Intergenic"
snps_df$Gene <- NA

# Check for gene overlaps
gene_hits <- findOverlaps(snp_gr, gene_gr)
snps_df$Annotation[queryHits(gene_hits)] <- "Genic"
snps_df$Gene[queryHits(gene_hits)] <- mcols(gene_gr)$proteinId[subjectHits(gene_hits)]

# Check for promotor overlaps (only if not already genic)
promotor_hits <- findOverlaps(snp_gr, promotor_gr)
not_genic <- snps_df$Annotation != "Genic"
snps_df$Annotation[queryHits(promotor_hits)[not_genic[queryHits(promotor_hits)]]] <- "promotor"
snps_df$Gene[queryHits(promotor_hits)[not_genic[queryHits(promotor_hits)]]] <- mcols(promotor_gr)$proteinId[subjectHits(promotor_hits)[not_genic[queryHits(promotor_hits)]]]

# call the df
snps_df
na.omit(unique(snps_df$Gene))

# Save annotated SNPs
write.csv(snps_df, "FST1_annotated_snps_with_genes.csv")


# now for top 1% of FST windows
# first for FST = 1 windows only

# gff and regions same as previous (so are gene_gr and promotor_gr)

# Read your FST regions
snps_df1 <- read.csv("snps_in_fst_top_1p_windows.csv")  

# Convert SNPs to GRanges
snp_gr1 <- GRanges(seqnames = snps_df1$CHROM,
                  ranges = IRanges(start = snps_df1$POS, end = snps_df1$POS),
                  SNP = snps_df1$SNP)

# Initialize annotation column
snps_df1$Annotation <- "Intergenic"
snps_df1$Gene <- NA

# Check for gene overlaps
# some scaffolds do not have annotated genes - this throws an error for findOverlaps()
# find the scaffolds only in both the snps and the gr sets
common_scaffolds <- intersect(seqlevels(snp_gr1), seqlevels(gene_gr))
snp_gr1 <- snp_gr1[seqnames(snp_gr1) %in% common_scaffolds]

# Restrict both GRanges objects to only the common scaffolds
snp_gr1 <- keepSeqlevels(snp_gr1, common_scaffolds, pruning.mode = "coarse")
gene_gr1 <- keepSeqlevels(gene_gr, common_scaffolds, pruning.mode = "coarse")

# Check for gene overlaps
gene_hits1 <- findOverlaps(snp_gr1, gene_gr1)
snps_df1$Annotation[queryHits(gene_hits1)] <- "Genic"
snps_df1$Gene[queryHits(gene_hits1)] <- mcols(gene_gr1)$proteinId[subjectHits(gene_hits1)]

# Check for promotor overlaps (only if not already genic)
# Restrict promotor GRanges objects to only the common scaffolds
promotor_gr1 <- keepSeqlevels(promotor_gr, common_scaffolds, pruning.mode = "coarse")

promotor_hits1 <- findOverlaps(snp_gr1, promotor_gr1)
not_genic1 <- snps_df1$Annotation != "Genic"
snps_df1$Annotation[queryHits(promotor_hits1)[not_genic1[queryHits(promotor_hits1)]]] <- "promotor"
snps_df1$Gene[queryHits(promotor_hits1)[not_genic1[queryHits(promotor_hits1)]]] <- mcols(promotor_gr1)$proteinId[subjectHits(promotor_hits1)[not_genic1[queryHits(promotor_hits1)]]]

# call the df
snps_df1
unique(snps_df1$Gene)

# Save annotated SNPs
write.csv(snps_df1, "FST_top_1p_annotated_snps_with_genes.csv")


# now for top 5% of FST windows
# first for FST = 1 windows only

# gff and regions same as previous (so are gene_gr and promotor_gr)

# Read your FST regions
snps_df5 <- read.csv("snps_in_fst_top_5p_windows.csv")  

# Convert SNPs to GRanges
snp_gr5 <- GRanges(seqnames = snps_df5$CHROM,
                   ranges = IRanges(start = snps_df5$POS, end = snps_df5$POS),
                   SNP = snps_df5$SNP)

# Initialize annotation column
snps_df5$Annotation <- "Intergenic"
snps_df5$Gene <- NA

# Check for gene overlaps
# some scaffolds do not have annotated genes - this throws an error for findOverlaps()
# find the scaffolds only in both the snps and the gr sets
common_scaffolds <- intersect(seqlevels(snp_gr5), seqlevels(gene_gr))
snp_gr5 <- snp_gr5[seqnames(snp_gr5) %in% common_scaffolds]

# Restrict both GRanges objects to only the common scaffolds
snp_gr5 <- keepSeqlevels(snp_gr5, common_scaffolds, pruning.mode = "coarse")
gene_gr5 <- keepSeqlevels(gene_gr, common_scaffolds, pruning.mode = "coarse")

# Check for gene overlaps
gene_hits5 <- findOverlaps(snp_gr5, gene_gr5)
snps_df5$Annotation[queryHits(gene_hits5)] <- "Genic"
snps_df5$Gene[queryHits(gene_hits5)] <- mcols(gene_gr5)$proteinId[subjectHits(gene_hits5)]

# Check for promotor overlaps (only if not already genic)
# Restrict promotor GRanges objects to only the common scaffolds
promotor_gr5 <- keepSeqlevels(promotor_gr, common_scaffolds, pruning.mode = "coarse")

promotor_hits5 <- findOverlaps(snp_gr5, promotor_gr5)
not_genic5 <- snps_df5$Annotation != "Genic"
snps_df5$Annotation[queryHits(promotor_hits5)[not_genic5[queryHits(promotor_hits5)]]] <- "Promotor"
snps_df5$Gene[queryHits(promotor_hits5)[not_genic5[queryHits(promotor_hits5)]]] <- mcols(promotor_gr5)$proteinId[subjectHits(promotor_hits5)[not_genic5[queryHits(promotor_hits5)]]]

# call the df
snps_df5
unique(snps_df5$Gene)

# Save annotated SNPs
write.csv(snps_df5, "FST_top_5p_annotated_snps_with_genes.csv")