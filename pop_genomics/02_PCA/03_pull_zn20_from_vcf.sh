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

singularity exec https://depot.galaxyproject.org/singularity/bcftools:1.9--h68d8f2e_9 bcftools view -S zn20.txt -o stomentosus_zn20_biallelic.vcf stomentosus_filtered_biallelic_only.vcf
