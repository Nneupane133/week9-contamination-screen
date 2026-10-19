#!/bin/bash
#PBS -N DNA_WGS_download
#PBS -q large
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o /home/aubnxn001/week9/logs/

# Load modules before enabling strict error handling.
set +e
set +u
set +E
set +o pipefail

source /apps/profiles/modules_asax.sh.dyn
module load sra/3.0.0

set -Eeuo pipefail
trap 'rc=$?; echo "ERROR: line $LINENO; exit $rc; command: $BASH_COMMAND" >&2; exit "$rc"' ERR

THREADS=8
ACCESSION="SRR39197043"
OUTDIR="/home/aubnxn001/week9/data/DNA_WGS"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Started: $(date)"

# 1. Download the .sra file into $OUTDIR/$ACCESSION/
prefetch "$ACCESSION" --output-directory "$OUTDIR" --max-size 50G --progress

# 2. Check the download is intact
vdb-validate "$OUTDIR/$ACCESSION"

# 3. Convert to paired FASTQ (_1 and _2), using $OUTDIR for temp space
fasterq-dump "$OUTDIR/$ACCESSION" \
    --split-files \
    --threads "$THREADS" \
    --temp "$OUTDIR" \
    --outdir "$OUTDIR" \
    --progress

# 4. Check both FASTQs were written
for f in "$OUTDIR/${ACCESSION}_1.fastq" "$OUTDIR/${ACCESSION}_2.fastq"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty $f" >&2; exit 1; }
done

# Compression and removal of the .sra file are done by compress_DNA_WGS.sh,
# so the .sra is kept here in case the FASTQs need to be regenerated.
ls -lh "$OUTDIR"
echo "Finished: $(date)"
echo "Next: qsub /home/aubnxn001/week9/scripts/compress_DNA_WGS.sh"
