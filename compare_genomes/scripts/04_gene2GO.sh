#!/bin/bash
#SBATCH --job-name=make_gene2go
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

set -euo pipefail

# ----------------------------- Paths (edit as needed) -----------------------------
# idmap: proteinId<TAB>jgi|... full header id used in Orthogroups.tsv
A_MAP="../idmap/isolateA.idmap.tsv"
B_MAP="../idmap/isolateB.idmap.tsv"

# GO tables (columns: proteinId go_term_id go_term_type go_name go_acc)
A_GO="../annotations/isolateA.GO.tsv"
B_GO="../annotations/isolateB.GO.tsv"

# Outputs
OUTDIR="../annotations"         # gene2go outputs
LOGDIR="../logs"                # coverage logs
MAKE_GLOBAL=1                   # set to 0 if you do not want global gene2go

mkdir -p "${LOGDIR}"

# ----------------------------- Functions -----------------------------------------

# Make one gene2go file from idmap + GO
# Args: idmap_path go_tsv out_gene2go iso_name
make_gene2go() {
  local idmap="$1"
  local gofile="$2"
  local out="$3"
  local iso="$4"

  if [[ ! -s "$idmap" ]]; then
    echo "ERROR: idmap missing or empty: $idmap" >&2; exit 1
  fi
  if [[ ! -s "$gofile" ]]; then
    echo "ERROR: GO file missing or empty: $gofile" >&2; exit 1
  fi

  # Join idmap (proteinId -> OrthoFinder ID) with GO rows (use go_acc column),
  # aggregate per gene, deduplicate GO terms, and output: gene<TAB>GO1,GO2,...
  awk 'BEGIN{FS=OFS="\t"}
       FNR==NR {map[$1]=$2; next}         # idmap
       NR==1 {next}                       # skip GO header
       {
         gid = map[$1];                   # $1 = proteinId in GO file
         go  = $5;                        # $5 = go_acc (GO:xxxx)
         if (gid != "" && go ~ /^GO:[0-9]+$/) {
           print gid, go
         }
       }' "$idmap" "$gofile" \
  | awk 'BEGIN{FS=OFS="\t"}
         {
           g=$1; term=$2;
           key=g "_" term;
           if (!(key in seen)) {
             if (list[g]=="") list[g]=term;
             else list[g]=list[g]","term;
             seen[key]=1;
           }
         }
         END{
           for (g in list) print g, list[g]
         }' \
  | sort -t$'\t' -k1,1 > "$out"

  # ----------------------------- Cleaning step (added) -----------------------------
  # After writing $out, sanitize to enforce exactly two columns and pure header ID:
  awk -F'\t' 'NF>=2{
    split($1, a, /[[:space:]]+/); gene=a[1];
    go=$2; gsub(/^[[:space:]]+|[[:space:]]+$/, "", gene); gsub(/^[[:space:]]+|[[:space:]]+$/, "", go);
    if (gene!="" && go!="") print gene "\t" go
  }' "$out" | sort -t$'\t' -k1,1 -u > "${out}.clean"
  mv "${out}.clean" "$out"
  # -------------------------------------------------------------------------------

  # Per-isolate coverage log
  local total_genes annotated_genes unique_go_terms total_go_rows pct
  total_genes=$(cut -f2 "$idmap" | sort -u | wc -l | awk '{print $1}')
  annotated_genes=$(cut -f1 "$out" | sort -u | wc -l | awk '{print $1}')
  unique_go_terms=$(cut -f2 "$out" | tr ',' '\n' | sort -u | wc -l | awk '{print $1}')
  total_go_rows=$(($(wc -l < "$gofile") - 1))
  pct=$(python3 - <<PY
t=$total_genes; a=$annotated_genes
print(f"{(a/t)*100:.2f}" if t>0 else "0.00")
PY
  )
  {
    echo "isolate: $iso"
    echo "idmap: $idmap"
    echo "gofile: $gofile"
    echo "gene2go: $out"
    echo "total_genes_in_fasta: $total_genes"
    echo "annotated_genes_with_GO: $annotated_genes"
    echo "unique_GO_terms: $unique_go_terms"
    echo "GO_rows(excl_header): $total_go_rows"
    echo "pct_annotated: $pct%"
  } > "${LOGDIR}/${iso}.gene2go.coverage.txt"
}

# ----------------------------- Run per isolate ----------------------------------

# A
make_gene2go "$A_MAP" "$A_GO" "${OUTDIR}/isolateA.gene2go.tsv" "isolateA"
# B
make_gene2go "$B_MAP" "$B_GO" "${OUTDIR}/isolateB.gene2go.tsv" "isolateB"

# ----------------------------- Optional global gene2go --------------------------
if [[ "$MAKE_GLOBAL" -eq 1 ]]; then
  cat "${OUTDIR}/isolateA.gene2go.tsv" \
      "${OUTDIR}/isolateB.gene2go.tsv" \
    | sort -u > "${OUTDIR}/global.gene2go.tsv"
fi

echo "Done. Outputs:"
echo " - ${OUTDIR}/isolateA.gene2go.tsv"
echo " - ${OUTDIR}/isolateB.gene2go.tsv"
[[ "$MAKE_GLOBAL" -eq 1 ]] && echo " - ${OUTDIR}/global.gene2go.tsv"
echo " - ${LOGDIR}/isolateA.gene2go.coverage.txt"
echo " - ${LOGDIR}/isolateB.gene2go.coverage.txt"

