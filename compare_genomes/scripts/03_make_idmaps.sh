#!/bin/bash
#SBATCH --job-name=orthofinder_3
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

A_PROT="../data/isolateA/proteins.faa"
B_PROT="../data/isolateB/proteins.faa"

OUTDIR="../idmap/"
mkdir -p "$OUTDIR"

# build ID maps for each genome
# gives proteinID from the full JGI name

grep "^>" "$A_PROT" | sed 's/>//' | awk -F"|" '{print $3"\t"$0}' > "$OUTDIR"/isolateA.idmap.tsv
grep "^>" "$B_PROT" | sed 's/>//' | awk -F"|" '{print $3"\t"$0}' > "$OUTDIR"/isolateB.idmap.tsv

echo "Done. Outputs:"
echo " - ${OUTDIR}/isolateA.idmap.tsv"
echo " - ${OUTDIR}/isolateB.idmap.tsv"