#!/bin/bash
# ---------------------------------------------------------------------------
# STEP 3a: fetch reference genomes for every candidate raised in steps 1-2,
#          plus a fixed panel of DECOYS.
#
# WHY DECOYS ARE NOT OPTIONAL:
#   If you map reads to a Chlamydia genome on its own, every conserved rRNA /
#   ribosomal-protein read in the library -- chicken, gut flora, anything --
#   will pile onto those conserved loci and manufacture a convincing-looking
#   "detection". The fix is COMPETITIVE mapping: put the host, the candidates
#   and a panel of likely alternative explanations into ONE reference, and let
#   each read go to its single best match. A read only counts for a candidate
#   if that candidate explains it better than the host or any decoy does.
#
#   The decoy panel is chosen to cover the alternative source categories we are
#   required to adjudicate between:
#     Homo sapiens      -> "human" (operator / index hopping from a human lane)
#     PhiX174           -> "lab organism" (Illumina spike-in)
#     E. coli K-12      -> "lab organism" / cloning host
#     S. cerevisiae     -> "lab organism" / common reagent eukaryote
#     UniVec_Core       -> "reagent" (cloning and adapter vector sequence)
#   Without these in the reference, reads from any of them would be forced onto
#   whichever candidate happened to look closest, and we would misattribute the
#   source even when the detection itself is real.
#
# Manifest: metadata/candidates.tsv, 4 tab-separated columns
#   name<TAB>role<TAB>taxid<TAB>url
#   role = host | candidate | decoy
# ---------------------------------------------------------------------------
set -euo pipefail
WK=/home/aubnxn001/week9
MAN=$WK/metadata/candidates.tsv
REF=$WK/references/candidates
mkdir -p "$REF"

[ -f "$MAN" ] || { echo "ERROR: manifest $MAN not found" >&2; exit 1; }

while IFS=$'\t' read -r name role taxid url; do
  case "$name" in ''|'#'*) continue;; esac
  dest="$REF/${name}.fna.gz"
  if [ -s "$dest" ]; then echo "  have $name"; continue; fi
  echo "[$(date +%T)] fetching $name ($role)"
  curl -fL --no-progress-meter -o "$dest" "$url" || { echo "  FAILED $name"; rm -f "$dest"; continue; }
  # UniVec ships as plain fasta; normalise everything to .fna.gz
  if ! gzip -t "$dest" 2>/dev/null; then mv "$dest" "${dest%.gz}"; gzip -f "${dest%.gz}"; fi
done < "$MAN"
echo "[$(date +%T)] candidate references:"; ls -la "$REF"
