#!/bin/bash
#SBATCH --job-name=make_sets
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan
# =============================================================================
# 05_build_sets.sh
# Build per-isolate gene sets (core/singleton/accessory) by joining:
#   - Orthogroups.tsv   (gene IDs per isolate)
#   - orthogroups_categories.tsv (presence flags + Category)
#
# Edit the two path variables below to match your file locations.
# =============================================================================

set -euo pipefail

# ----------------------------- Paths (EDIT THESE) -----------------------------
# Path to OrthoFinder's Orthogroups.tsv (tab-delimited; columns: Orthogroup, isolateA, isolateB, isolateC)
ORTHOGROUPS_TSV="../results/orthofinder/input/OrthoFinder/Results_Dec15/Orthogroups/Orthogroups.tsv"

# Path to your categories table (tab-delimited; columns: Orthogroup, isolateA, isolateB, isolateC, Category)
CATEGORIES_TSV="../results/core_accessory_singleton/orthogroup_categories.tsv"

# ----------------------------- Output dirs -----------------------------------
mkdir -p ../results/sets/core
mkdir -p ../results/sets/singleton
mkdir -p ../results/sets/accessory


# ----------------------------- Sanity checks ---------------------------------
if [[ ! -s "$ORTHOGROUPS_TSV" ]]; then
  echo "ERROR: Orthogroups.tsv not found or empty: $ORTHOGROUPS_TSV" >&2
  exit 1
fi
if [[ ! -s "$CATEGORIES_TSV" ]]; then
  echo "ERROR: orthogroups_categories.tsv not found or empty: $CATEGORIES_TSV" >&2
  echo "TIP: Set CATEGORIES_TSV to the correct location of your categories file." >&2
  exit 1
fi

# ----------------------------- Normalize CRLF (Windows line endings) ----------
# Only run conversion if files exist
if grep -q $'\r' "$ORTHOGROUPS_TSV"; then
  echo "INFO: Converting CRLF to LF in $ORTHOGROUPS_TSV"
  perl -pi -e 's/\r$//' "$ORTHOGROUPS_TSV"
fi
if grep -q $'\r' "$CATEGORIES_TSV"; then
  echo "INFO: Converting CRLF to LF in $CATEGORIES_TSV"
  perl -pi -e 's/\r$//' "$CATEGORIES_TSV"
fi

# ----------------------------- Build sets ------------------------------------
awk '
BEGIN { FS = "\t"; OFS = "\t" }

# ---------- Pass 1: read Orthogroups.tsv and store per-isolate gene strings ----------
FNR == NR {
  if (NR == 1) {
    # Header validation
    idxA = idxB = idxC = 0
    for (i = 1; i <= NF; i++) {
      if ($i == "isolateA") idxA = i
      else if ($i == "isolateB") idxB = i
    }
    if (!idxA || !idxB) {
      print "ERROR: Missing isolateA/isolateB columns in Orthogroups.tsv header." > "/dev/stderr"
      exit 1
    }
    next
  }
  og = $1
  a = $idxA; b = $idxB

  # Trim whitespace
  gsub(/^[ \t]+|[ \t]+$/, "", a)
  gsub(/^[ \t]+|[ \t]+$/, "", b)

  # Store raw strings (may be empty)
  Araw[og] = a; Braw[og] = b
  next
}

# ---------- Pass 2: read orthogroups_categories.tsv and emit set files ----------
FNR == 1 {
  # Header mapping
  ogIdx = aFlagIdx = bFlagIdx = catIdx = 0
  for (i = 1; i <= NF; i++) {
    if     ($i == "Orthogroup") ogIdx    = i
    else if($i == "isolateA")   aFlagIdx = i
    else if($i == "isolateB")   bFlagIdx = i
    else if($i == "Category")   catIdx   = i
  }
  if (!ogIdx || !aFlagIdx || !bFlagIdx || !catIdx) {
    print "ERROR: Missing Orthogroup/isolateA/isolateB/Category columns in categories header." > "/dev/stderr"
    exit 1
  }
  coreOG = singOG = 0
  next
}

{
  og    = $ogIdx
  aflag = $aFlagIdx
  bflag = $bFlagIdx
  cat   = $catIdx

  # Normalize flags and category text
  gsub(/^[ \t]+|[ \t]+$/, "", aflag)
  gsub(/^[ \t]+|[ \t]+$/, "", bflag)
  gsub(/^[ \t]+|[ \t]+$/, "", cat)

  # Pull gene strings from pass 1
  A = Araw[og]; B = Braw[og]; C = Craw[og]
  if (!(og in Araw) && !(og in Braw) && !(og in Craw)) {
    print "WARN: Orthogroup " og " present in categories but missing in Orthogroups.tsv" > "/dev/stderr"
    next
  }

  # Split lists by comma + optional space
  nA = split(A, arrA, /,[ ]*/)
  nB = split(B, arrB, /,[ ]*/)

  # Category handling (case-insensitive without relying on tolower)
  catL = cat
  if      (cat == "CORE"      || cat == "core")      catL = "core"
  else if (cat == "SINGLETON" || cat == "singleton") catL = "singleton"

  if (catL == "core") {
    coreOG++
    # per-isolate core sets
    for (i = 1; i <= nA; i++) { g = arrA[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/core/isolateA.core.genes.txt" }
    for (i = 1; i <= nB; i++) { g = arrB[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/core/isolateB.core.genes.txt" }

    # global core union
    for (i = 1; i <= nA; i++) { g = arrA[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/core/global.core.genes.txt" }
    for (i = 1; i <= nB; i++) { g = arrB[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/core/global.core.genes.txt" }
  }
  else if (catL == "singleton") {
    singOG++
    # Only the isolate with flag==1 contributes genes
    if (aflag == "1") for (i = 1; i <= nA; i++) { g = arrA[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/singleton/isolateA.singleton.genes.txt" }
    if (bflag == "1") for (i = 1; i <= nB; i++) { g = arrB[i]; gsub(/^[ \t]+|[ \t]+$/, "", g); if (g != "") print g >> "../results/sets/singleton/isolateB.singleton.genes.txt" }
  }
  else {
    print "WARN: Unrecognized Category for " og ": " cat > "/dev/stderr"
  }
}
END {
  print "INFO: CORE OGs: "      coreOG > "/dev/stderr"
  print "INFO: SINGLETON OGs: " singOG > "/dev/stderr"
}
' "$ORTHOGROUPS_TSV" "$CATEGORIES_TSV"

# ----------------------------- Deduplicate & sort --------------------------------
# Use find to avoid glob errors when no files exist yet
find ../results/sets -type f -name "*.genes.txt" -print0 | xargs -0 -r -I{} bash -c 'sort -u "{}" -o "{}"'

# ----------------------------- Log counts ---------------------------------------
find ../results/sets -type f -name "*.genes.txt" -print0 \
  | xargs -0 -r wc -l \
  | tee ../logs/set_counts.txt
