#!/bin/bash
#SBATCH --job-name=coverage_sets
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

set -euo pipefail

# annotation gene2go files (made in 04_gene2GO.sh)
GENE2GO_DIR="../annotations"
# set files (from 05_make_sets.sh)
SETS_DIR="../results/sets"
# output dir for log files
LOG_DIR="../logs"

# function to summarize isoaltes
summarize_isolate() {
  local iso="$1"
  local gene2go="${GENE2GO_DIR}/${iso}.gene2go.tsv"
  local out="${LOG_DIR}/coverage_summary.tsv"

  if [[ ! -s "$gene2go" ]]; then
    echo "WARN: missing gene2go for ${iso}: ${gene2go}" >&2
    # still write zeros
    for cat in core singleton accessory; do
      printf "%s\t%s.%s\t%d\t%d\t%d\t%s\n" "$iso" "$iso" "$cat" 0 0 0 "0.00" >> "$out"
    done
    return 0
  fi

  local ann_list
  ann_list="$(mktemp)"
  cut -f1 "$gene2go" | sort -u > "$ann_list"

  for cat in core singleton accessory; do
    local setfile="${SETS_DIR}/${cat}/${iso}.${cat}.genes.txt"
    local setname="${iso}.${cat}"
    if [[ ! -e "$setfile" ]]; then
      echo "WARN: setfile not found: $setfile" >&2
      printf "%s\t%s\t%d\t%d\t%d\t%s\n" "$iso" "$setname" 0 0 0 "0.00" >> "$out"
      continue
    fi

    local total annotated unannotated pct
    total=$(sort -u "$setfile" | wc -l | awk '{print $1}')
    if [[ "$total" -eq 0 ]]; then
      printf "%s\t%s\t%d\t%d\t%d\t%s\n" "$iso" "$setname" 0 0 0 "0.00" >> "$out"
      continue
    fi
    annotated=$(comm -12 <(sort -u "$setfile") "$ann_list" | wc -l | awk '{print $1}')
    unannotated=$(( total - annotated ))
    pct=$(python3 - <<PY
t=$total; a=$annotated
print(f"{(a/t)*100:.2f}" if t>0 else "0.00")
PY
    )
    printf "%s\t%s\t%d\t%d\t%d\t%s\n" "$iso" "$setname" "$total" "$annotated" "$unannotated" "$pct" >> "$out"
  done

  rm -f "$ann_list"
}

# echo for slurm output - where files are going
echo -e "isolate\tset\ttotal_genes\tannotated_genes\tunannotated_genes\tpct_annotated" > "${LOG_DIR}/coverage_summary.tsv"

# run summary function for each isolate
summarize_isolate "isolateA"
summarize_isolate "isolateB"

# global core summary
if [[ -s "${GENE2GO_DIR}/global.gene2go.tsv" && -s "${SETS_DIR}/core/global.core.genes.txt" ]]; then
  ann_list="$(mktemp)"
  cut -f1 "${GENE2GO_DIR}/global.gene2go.tsv" | sort -u > "$ann_list"

  total=$(sort -u "${SETS_DIR}/core/global.core.genes.txt" | wc -l | awk '{print $1}')
  if [[ "$total" -eq 0 ]]; then
    printf "%s\t%s\t%d\t%d\t%d\t%s\n" "GLOBAL" "global.core" 0 0 0 "0.00" >> "${LOG_DIR}/coverage_summary.tsv"
  else
    annotated=$(comm -12 <(sort -u "${SETS_DIR}/core/global.core.genes.txt") "$ann_list" | wc -l | awk '{print $1}')
    unannotated=$(( total - annotated ))
    pct=$(python3 - <<PY
t=$total; a=$annotated
print(f"{(a/t)*100:.2f}" if t>0 else "0.00")
PY
    )
    printf "%s\t%s\t%d\t%d\t%d\t%s\n" "GLOBAL" "global.core" "$total" "$annotated" "$unannotated" "$pct" >> "${LOG_DIR}/coverage_summary.tsv"
  fi
  rm -f "$ann_list"
fi

echo "Coverage summary written: ${LOG_DIR}/coverage_summary.tsv"
