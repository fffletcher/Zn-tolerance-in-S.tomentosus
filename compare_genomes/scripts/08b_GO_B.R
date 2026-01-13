
#!/usr/bin/env Rscript

# GO enrichment with topGO (BP, MF, CC) using Fisher's exact test (classic, elim, weight)
# Input:
#   - gene2go.tsv: gene<TAB>GO1,GO2,...
#   - foreground gene sets: one or more files with one gene ID per line (no header)
#   - background universe: optional file with one gene ID per line; defaults to genes in gene2go
# Output:
#   - TSVs per ontology and per gene set with enrichment results
#   - Optional PDFs with top terms

library(topGO)

setwd("~/Desktop/compGenome2/")

# ---------------------------- CONFIG (edit as needed) ----------------------------
# Example paths; replace with isolate-specific files or pass via a wrapper script


# generate files for accessory genes with ZnT names, for each group
# accessory A&B not C and A&C not B
keep_prefix <- function(infile, outfile, prefix_regex = "^jgi\\|Suito220_") {
  x <- readLines(infile, warn = FALSE)
  # If multiple IDs per line, split on commas/semicolons/whitespace
  x <- unlist(strsplit(x, "[,;[:space:]]+"))
  x <- trimws(x)
  x <- x[grepl(prefix_regex, x)]
  x <- sort(unique(x[x != ""]))
  writeLines(x, outfile)
}


GENE2GO_FILE <- "annotations/isolateB.gene2go.tsv"  # cleaned two-column file
FOREGROUND_FILES <- c(
  "results/sets/core/isolateB.core.genes.txt",
  "results/sets/singleton/isolateB.singleton.genes.txt"
)
BACKGROUND_FILE <- NA  # set to a path to override default background; e.g., "sets/isolateA_universe.txt"

OUTDIR <- "results/go_enrichment/isolateB"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE) # outputs go here
ALGORITHMS <- c("weight")  # topGO algorithms to run
ONTOLOGIES <- c("BP", "MF", "CC")             # biological process, molecular function, cellular component

TOP_N_TO_PLOT <- 15   # top terms to plot per set+ontology (PDF)
MIN_TERM_SIZE <- 5    # min size of GO term to keep (filters extremely small terms)
PVAL_CUTOFF <- 0.05   # report table rows with FDR <= this

# ---------------------------- Helpers -------------------------------------------

read_gene2go <- function(path) {
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("gene", "terms_csv")
  # Split comma-separated terms into a list
  terms_list <- strsplit(x$terms_csv, ",")
  names(terms_list) <- x$gene
  # Ensure GO IDs valid
  terms_list <- lapply(terms_list, function(v) {
    v <- trimws(v)
    v[grepl("^GO:[0-9]+$", v)]
  })
  # Drop genes with no valid terms
  terms_list <- terms_list[lengths(terms_list) > 0]
  return(terms_list)
}

read_gene_list <- function(path) {
  if (is.na(path)) return(character(0))
  v <- readLines(path, warn = FALSE)
  v <- unique(trimws(v))
  v <- v[v != ""]
  return(v)
}

build_universe <- function(gene2go_list, background_file = NA) {
  if (is.na(background_file)) {
    # Default: all genes present in gene2go mapping
    return(sort(unique(names(gene2go_list))))
  } else {
    user_bg <- read_gene_list(background_file)
    # Intersect with genes that have GO terms (topGO requirement)
    return(sort(intersect(user_bg, names(gene2go_list))))
  }
}

make_gene_selection <- function(universe, foreground) {
  # Named numeric vector 1/0 with names = universe genes
  sel <- rep(0, length(universe))
  names(sel) <- universe
  sel[names(sel) %in% foreground] <- 1
  return(sel)
}

run_topgo_for_set <- function(set_name, universe, geneSel, gene2go, ontology) {
  # Build topGOdata object
  obj <- new("topGOdata",
             ontology = ontology,
             allGenes = geneSel,
             geneSelectionFun = function(x) x == 1,
             annot = annFUN.gene2GO,
             gene2GO = gene2go,
             nodeSize = MIN_TERM_SIZE)
  
  results <- list()
  
  # Fisher's exact test across algorithms
  for (alg in ALGORITHMS) {
    test_stat <- runTest(obj, algorithm = alg, statistic = "fisher")
    tab <- GenTable(obj, 
                    p.value = test_stat,
                    topNodes = length(usedGO(obj)),
                    numChar = 100)
    # Add metadata columns
    tab$algorithm <- alg
    tab$set <- set_name
    tab$ontology <- ontology
    
    # Add raw p.values (extract from test_stat slot) + FDR
    # The p-values in GenTable can be rounded; grab exact values:
    pvals <- score(test_stat)  # named vector by GO term
    tab$raw_p <- pvals[tab$GO.ID]
    tab$raw_p[is.na(tab$raw_p)] <- tab$p.value[is.na(tab$raw_p)]  # fallback
    tab$FDR <- p.adjust(tab$raw_p, method = "BH")
    
    # Order by FDR, then raw p
    tab <- tab[order(tab$FDR, tab$raw_p, decreasing = FALSE), ]
    results[[alg]] <- tab
  }
  
  # Combine across algorithms
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  return(list(topGOobj = obj, table = out))
}


plot_top_terms <- function(res_table, set_name, ontology, outdir,
                           top_n = 20, fdr_plot = 0.25) {
  # Prefer weight terms that pass the FDR threshold; fallback to classic if none
  sel <- subset(res_table, algorithm == "weight" & FDR <= fdr_plot)
  if (nrow(sel) < 1) sel <- subset(res_table, algorithm == "classic")
  sel <- head(sel, top_n)
  if (nrow(sel) < 1) return(invisible(NULL))
  
  # Make sure output directory exists
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  pdf(file.path(outdir, sprintf("%s_%s_top_terms.pdf", set_name, ontology)),
      width = 7, height = 9)
  vals <- -log10(sel$FDR)
  names(vals) <- paste0(sel$GO.ID, " ", sel$Term)
  
  op <- par(mar = c(8, 12, 4, 2))
  barplot(vals,
          horiz = TRUE,
          las = 1,
          cex.names = 0.7,
          main = sprintf("%s (%s): topGO top terms", set_name, ontology),
          xlab = expression(-log[10]~FDR))
  par(op)
  dev.off()
}


# ---------------------------- Main ----------------------------------------------

gene2go <- read_gene2go(GENE2GO_FILE)
cat(sprintf("Loaded %d genes with GO from %s\n", length(gene2go), GENE2GO_FILE))

universe <- build_universe(gene2go, BACKGROUND_FILE)
cat(sprintf("Universe size: %d genes\n", length(universe)))

# Filter gene2go to universe only
gene2go <- gene2go[names(gene2go) %in% universe]

# Read all foreground gene sets
sets <- list()
for (f in FOREGROUND_FILES) {
  sname <- tools::file_path_sans_ext(basename(f))
  genes <- read_gene_list(f)
  # Keep only genes in universe
  genes <- intersect(genes, universe)
  sets[[sname]] <- sort(unique(genes))
  cat(sprintf("Set %s: %d genes (after universe intersect)\n", sname, length(sets[[sname]])))
}

# Iterate over sets and ontologies
# Run enrichment per set and ontology
for (sname in names(sets)) {
  sel <- make_gene_selection(universe, sets[[sname]])
  for (ont in ONTOLOGIES) {
    res <- run_topgo_for_set(set_name = sname, universe = universe, geneSel = sel,
                             gene2go = gene2go, ontology = ont)
    tab <- res$table
    
    # Write outputs
    out_full <- file.path(OUTDIR, sprintf("%s_%s_topGO_full.tsv", sname, ont))
    out_sig  <- file.path(OUTDIR, sprintf("%s_%s_topGO_sig.tsv", sname, ont))
    write.table(tab, out_full, sep = "\t", row.names = FALSE, quote = FALSE)
    
    sigtab <- subset(tab, FDR <= PVAL_CUTOFF)
    write.table(sigtab, out_sig, sep = "\t", row.names = FALSE, quote = FALSE)
    
    plot_top_terms(tab, sname, ont, OUTDIR, top_n = TOP_N_TO_PLOT, fdr_plot = 0.25)
    
    cat(sprintf("Saved: %s and %s (%d significant at FDR <= %.2f)\n",
                out_full, out_sig, nrow(sigtab), PVAL_CUTOFF))
  }
}
