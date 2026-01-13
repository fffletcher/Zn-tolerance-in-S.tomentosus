# Zn-tolerance-in-S.tomentosus
Data and scripts for analysis xxxxxxx

compare_genomes scripts - 
01_orthofinder.sh = use othofinder to generate orthogroups between genomes A and B
02_orthogroups.sh = orthogroup presence/absence per genome, and which genes are in each orthogroup
03_make_idmaps.sh = generate table that gives gene name from JGI naming conventions
04_gene2GO.sh = make gene2go file from annotations - to be used in topGo for GO enrichment analysis
05_make_sets.sh = generate lists of genes for each isolate that are either core or singleton genes
06_coverage_sets.sh = check genes in the gff file (generated from PacBio reads) have coverage in short read sequencing (Illumina reads)

- BUSCO: Get BUSCO scores for each genome (ZnT/St230 = BUSCO_A.sh and ZnS/St220 = BUSCO_B.sh), generate Venn diagram of orthogroup presence/absence (Venn.R), then perform GO analysis via TopGO (GO_A.R, GO_B.R).
