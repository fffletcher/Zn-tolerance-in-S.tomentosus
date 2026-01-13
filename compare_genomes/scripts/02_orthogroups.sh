#!/bin/bash
#SBATCH --job-name=orthogroups
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

set -euo pipefail

ORTHO_RES="../results/orthofinder/input/OrthoFinder/Results_Dec15"
OG_DIR="${ORTHO_RES}/Orthogroups"
OG_COUNTS="${OG_DIR}/Orthogroups.GeneCount.tsv"
OG_TSV="${OG_DIR}/Orthogroups.tsv"

# Output directory
OUTDIR="../results/core_accessory_singleton"
mkdir -p "$OUTDIR"

require() { [[ -s "$1" ]] || { echo "ERROR: Missing file: $1" >&2; exit 1; }; }
hasfile() { [[ -n "${1:-}" && -s "$1" ]]; }

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Validate required inputs
require "$OG_COUNTS"
require "$OG_TSV"

echo "[OG/Gene] Using:"
echo "  Orthogroups counts: $OG_COUNTS"
echo "  Orthogroups genes : $OG_TSV"
echo "  Output dir         : $OUTDIR"

# Detect isolate names from OG_COUNTS header
ISO_NAMES_FILE="${OUTDIR}/isolate_names.tsv"
awk 'BEGIN{FS="\t"; OFS="\t"} NR==1 { 
  # header expected as: Orthogroup <tab> iso1 <tab> iso2 <tab> iso3
  if (NF < 4) { print "ERROR: Expected at least 4 columns in header" > "/dev/stderr"; exit 1 }
  print $2, $3, $4 
}' "$OG_COUNTS" > "$ISO_NAMES_FILE"

ISO1_NAME=$(awk 'NR==1{print $1}' "$ISO_NAMES_FILE")
ISO2_NAME=$(awk 'NR==1{print $2}' "$ISO_NAMES_FILE")

echo "[OG/Gene] Detected isolates: $ISO1_NAME, $ISO2_NAME"

# ================
# Step 2: Build orthogroup-level presence (1/0) and categories
# ================
OG_PRES="${OUTDIR}/orthogroup_presence.tsv"
OG_CAT="${OUTDIR}/orthogroup_categories.tsv"
OG_SUM="${OUTDIR}/orthogroup_category_counts.tsv"

awk -v iso1="$ISO1_NAME" -v iso2="$ISO2_NAME" '
  BEGIN{FS="\t"; OFS="\t"; c1=0;c2=0;c3=0; ogcol=0}
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i ~ /^Orthogroup$/) ogcol=i
      if ($i == iso1) c1=i
      if ($i == iso2) c2=i
    }
    if (ogcol==0) ogcol=1
    if (c1==0 || c2==0) { print "ERROR: Could not match isolate names in header" > "/dev/stderr"; exit 1 }
    print "Orthogroup", iso1, iso2
    next
  }
  {
    a=$c1+0; b=$c2+0; c=$c3+0;
    pa=(a>0)?1:0; pb=(b>0)?1:0; pc=(c>0)?1:0;
    print $ogcol, pa, pb
  }
' "$OG_COUNTS" > "$OG_PRES"

awk '
  BEGIN{FS="\t"; OFS="\t"}
  NR==1 { print $0, "Category"; next }
  {
    a=$2+0; b=$3+0;
    pres=a+b;
    cat=(pres==2)?"CORE":(pres==1)?"SINGLETON":"ABSENT";
    print $0, cat;
  }
' "$OG_PRES" > "$OG_CAT"

awk '
  BEGIN{FS="\t"}
  NR>1 { cnt[$NF]++ }
  END{ for(k in cnt) print k"\t"cnt[k] }
' "$OG_CAT" | sort > "$OG_SUM"

echo "[OG/Gene] Wrote:"
echo "  Orthogroup presence   : $OG_PRES"
echo "  Orthogroup categories : $OG_CAT"
echo "  Category counts       : $OG_SUM"


# Build gene-level table (Orthogroup, Isolate, Gene)
GENE_TABLE="${OUTDIR}/gene_table.tsv"

awk -v iso1="$ISO1_NAME" -v iso2="$ISO2_NAME" '
  BEGIN{FS="\t"; OFS="\t"}
  NR==1 { next }  # skip header line in Orthogroups.tsv
  {
    og=$1;
    # Column order in Orthogroups.tsv mirrors GeneCount.tsv header order
    # Split by comma+space
    n1=split($2,a1,", "); for(i=1;i<=n1;i++) if(a1[i]!="") print og, iso1, a1[i], 1;
    n2=split($3,a2,", "); for(i=1;i<=n2;i++) if(a2[i]!="") print og, iso2, a2[i], 1;
  }
' "$OG_TSV" > "$GENE_TABLE.raw"

# Add header
printf "Orthogroup\tIsolate\tGene\tPresent\n" > "$GENE_TABLE"
cat "$GENE_TABLE.raw" >> "$GENE_TABLE"
rm -f "$GENE_TABLE.raw"

echo "[OG/Gene] Wrote gene-level table: $GENE_TABLE"

echo "[OG/Gene] Done."
echo "[OG/Gene] Preview (categories):"
head -n 10 "$OG_CAT"
