#!/bin/bash
#SBATCH --job-name=BUSCO
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

set -euo pipefail

# Must be connected to internet to download the lineage

CONTAINER="http://depot.galaxyproject.org/singularity/busco:6.0.0--pyhdfd78af_0"
LINEAGE="basidiomycota_odb10"
CPUS=16
MODE=genome

# Find the genome FASTA file in the directory
INPUT="../data/isolateB/genome.fa"

if [[ -z "$INPUT" ]]; then
  echo "ERROR: No genome.fa file found in current directory."
  exit 1
fi

# set and make output directory
OUTDIR="../results/busco/isolateB"

mkdir -p "${OUTDIR}"

echo "Running BUSCO on ${INPUT} (mode=${MODE}, lineage=${LINEAGE})"

singularity exec \
  --bind "$(realpath ../data)":/data,"$(realpath ../results)":/results \
  "${CONTAINER}" \
  busco \
    --in "/data/isolateB/genome.fa" \
    --mode "${MODE}" \
    --lineage "${LINEAGE}" \
    --out "isolateA" \
    --out_path "/results/busco/isolateB" \
    --cpu "${CPUS}"

echo "BUSCO finished."
