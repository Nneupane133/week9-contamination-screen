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
# Manifest: metadata/candidates.tsv, tab separated
#   name  role  taxid  source  id  restrict
#   source = assembly (NCBI datasets) | nuccore (efetch) | url | local
# ---------------------------------------------------------------------------
set -euo pipefail
WK=/home/aubnxn001/week9
MAN=$WK/metadata/candidates.tsv
REF=$WK/references/candidates
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$REF"
set +e; source /apps/profiles/modules_asax.sh.dyn >/dev/null 2>&1
module load ncbi-datasets-cli seqkit >/dev/null 2>&1; set -e
EU=https://eutils.ncbi.nlm.nih.gov/entrez/eutils

[ -f "$MAN" ] || { echo "ERROR: manifest $MAN not found" >&2; exit 1; }

while IFS=$'\t' read -r name role taxid source id restrict; do
  case "$name" in ''|'#'*) continue;; esac
  [ "$source" = local ] && continue              # hosts are fetched by 00_fetch_host_refs.sh
  dest="$REF/${name}.fna.gz"
  if [ -s "$dest" ]; then echo "  have $name"; continue; fi
  echo "[$(date +%T)] fetching $name  ($role, $source $id)"
  case "$source" in
    assembly)
      # NCBI datasets gives a checksum-verified assembly and resolves the
      # versioned FTP path for us, which hand-built URLs get wrong constantly.
      ( cd "$TMP" && rm -rf ds && \
        datasets download genome accession "$id" --include genome --filename ds.zip \
          >/dev/null 2>&1 && unzip -q -o ds.zip -d ds ) || { echo "  FAILED $name"; continue; }
      cat "$TMP"/ds/ncbi_dataset/data/"$id"/*.fna 2>/dev/null | gzip > "$dest" || true
      ;;
    nuccore)
      # One or more comma-separated accessions (segmented viruses need several).
      : > "$TMP/$name.fa"
      ok=1
      for acc in ${id//,/ }; do
        curl -fsSg "$EU/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
          >> "$TMP/$name.fa" || ok=0
        sleep 0.4
      done
      [ "$ok" = 1 ] && [ -s "$TMP/$name.fa" ] && gzip -c "$TMP/$name.fa" > "$dest"
      ;;
    url)
      curl -fsSL -o "$TMP/$name.raw" "$id" || { echo "  FAILED $name"; continue; }
      if gzip -t "$TMP/$name.raw" 2>/dev/null; then cp "$TMP/$name.raw" "$dest"
      else gzip -c "$TMP/$name.raw" > "$dest"; fi
      ;;
  esac
  if [ -s "$dest" ]; then echo "  ok  $(zcat "$dest" | grep -c '^>') sequences"
  else echo "  FAILED $name"; rm -f "$dest"; fi
done < "$MAN"

echo; echo "[$(date +%T)] candidate references:"
for f in "$REF"/*.fna.gz; do printf "  %-34s %s\n" "$(basename "$f" .fna.gz)" "$(du -h "$f"|cut -f1)"; done
