#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

set -euo pipefail

module purge
module load singularity

A_PROT="../data/isolateA/proteins.faa"
B_PROT="../data/isolateB/proteins.faa"

INDIR="../results/orthofinder/input"
OUTDIR="../results/orthofinder"

mkdir -p "$INDIR"
mkdir -p "$OUTDIR"

# Copy the .fa files to the input directory
cp "$A_PROT" "$INDIR/isolateA.faa"
cp "$B_PROT" "$INDIR/isolateB.faa"


IMG_ORTHOFINDER="https://depot.galaxyproject.org/singularity/orthofinder:2.5.4--hdfd78af_0"


# Run Orthofinder (-t threads for BLAST/DIAMOND, -a threads for OrthoFinder steps)
singularity exec "$IMG_ORTHOFINDER" \
  orthofinder -f "$INDIR" -t "${SLURM_CPUS_PER_TASK}" -a "${SLURM_CPUS_PER_TASK}"
