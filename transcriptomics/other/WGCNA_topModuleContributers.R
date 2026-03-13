# pulling genes from WGCNA modules based on 
# significance for trait and module membership

module <- "green" 
column <- match(module, modNames)
moduleGenes <- moduleLabels==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceZn[moduleGenes,1]), 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance in Zn Treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Logical vector for genes in the module
moduleGenes <- moduleLabels == module

# Extract MM and GS values for the module
MM_values <- abs(geneModuleMembership[moduleGenes, column])
GS_values <- abs(geneTraitSignificanceZn[moduleGenes, 1])

# Calculate the 99th percentile thresholds
MM_cutoff <- quantile(MM_values, 0.90)
GS_cutoff <- quantile(GS_values, 0.90)

# Filter genes above both thresholds
topGenesLogical <- MM_values > MM_cutoff & GS_values > GS_cutoff

# Get gene names (assuming rownames are gene IDs)
topGenes <- rownames(geneModuleMembership)[moduleGenes][topGenesLogical]

# Optionally rank by combined score
combinedScore <- MM_values[topGenesLogical] * GS_values[topGenesLogical]
topGenesRanked <- topGenes[order(combinedScore, decreasing = TRUE)]
topGenesRanked
