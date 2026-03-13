#!/bin/bash
#SBATCH --time=0-6:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

module purge
module load singularity

# Using PLINK 1.9 
# Convert VCF file to PLINK’s .bed, .bim, and .fam formats:

singularity exec https://depot.galaxyproject.org/singularity/plink:1.90b7.7--h18e278d_1 plink --vcf stomentosus_filtered_zn20.vcf --allow-extra-chr --make-bed --out plink_zn20_output

# Perform PCA (generate eigenvec and eigenval files)

singularity exec https://depot.galaxyproject.org/singularity/plink:1.90b7.7--h18e278d_1 plink --bfile plink_zn20_output --allow-extra-chr --pca --out pca_zn20_results
