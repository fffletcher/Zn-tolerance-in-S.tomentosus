#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

module purge
module load singularity

singularity exec https://depot.galaxyproject.org/singularity/vcftools:0.1.17--pl5321h077b44d_0 vcftools --vcf stomentosus_filtered_biallelic_only.vcf --fst-window-size 5000 --weir-fst-pop top10.txt --weir-fst-pop bottom10.txt --out top_bottom.fst


