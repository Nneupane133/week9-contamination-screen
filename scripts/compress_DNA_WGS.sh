#!/bin/bash
#PBS -N WGS_gzip
#PBS -q large
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: line $LINENO; exit $rc; command: $BASH_COMMAND" >&2; exit "$rc"' ERR

THREADS=8
ACCESSION="SRR39197043"
OUTDIR="/home/aubnxn001/week9/data/DNA_WGS"

cd "$OUTDIR"
echo "Started: $(date) on $(hostname)"

# Stop before deleting anything if the uncompressed FASTQs are not here
# (e.g. ACCESSION is wrong or this sample was already compressed).
for f in "${ACCESSION}"_1.fastq "${ACCESSION}"_2.fastq; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty $OUTDIR/$f" >&2; exit 1; }
done

# Remove partial .gz files left by killed runs; the .fastq originals are intact.
rm -f "${ACCESSION}"_1.fastq.gz "${ACCESSION}"_2.fastq.gz

if command -v pigz >/dev/null 2>&1; then
    pigz -p "$THREADS" "${ACCESSION}"_1.fastq "${ACCESSION}"_2.fastq
else
    gzip "${ACCESSION}"_1.fastq "${ACCESSION}"_2.fastq
fi

# Remove the .sra file once FASTQs are verified
gzip -t "${ACCESSION}"_1.fastq.gz "${ACCESSION}"_2.fastq.gz
rm -rf "$ACCESSION"

ls -lh
echo "Finished: $(date)"
