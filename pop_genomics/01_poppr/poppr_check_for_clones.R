#load libraries
library(poppr)
library(adegenet)
library(vcfR)

#set seed for reproducibility
set.seed(1916)

#read in vcf
vcf_me<- read.vcfR("~/Downloads/stomentosus_filtered.vcf.gz") # filtered vcf from RepAdapt - contains multiallelic SNPS

# multiallelic snps filtered out in below step
#convert vcf to genlight object
genlight <- vcfR2genlight(vcf_me)
rm(vcf_me)

# 4043408 SNPs total
# 96498 multiallelic SNPs filtered out
biallelicSNPs <- 4043408 - 96498
biallelicSNPs

# use mlg function in poppr to identify clones
# Investigate the number of multilocus genotypes.
amlg <- mlg(genlight)
amlg 

# number of multilocus genotypes = 48
# same number as samples = no clones