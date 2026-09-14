#!/bin/bash
# rename_demux.sh
# Renames demultiplexed FASTQ files from barcode-pair names to sample IDs.
#
# Usage:
#   bash rename_demux.sh <sample_map.txt> [fastq_dir]
#
# Arguments:
#   sample_map.txt  : two-column file (tab or space separated):
#                     column 1 = barcode pair (e.g. ACGTTACC+TTGTCGGT)
#                     column 2 = sample name/ID
#   fastq_dir       : directory containing the .fastq.gz files (default: current dir)
#
# Expected input filename format : sample_<BARCODE>_R1.fastq.gz / _R2.fastq.gz
# Output filename format         : sample_<SAMPLEID>_R1.fastq.gz / _R2.fastq.gz

# ── argument handling ────────────────────────────────────────────────────────
SAMPLE_MAP="$1"
FASTQ_DIR="${2:-.}"

if [[ -z "$SAMPLE_MAP" ]]; then
    echo "Usage: $0 <sample_map.txt> [fastq_dir]"
    exit 1
fi

if [[ ! -f "$SAMPLE_MAP" ]]; then
    echo "ERROR: sample map not found: $SAMPLE_MAP"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

# ── counters ─────────────────────────────────────────────────────────────────
renamed=0
missing=0
skipped=0

echo "=================================================="
echo " Renaming demultiplexed FASTQ files"
echo " Map   : $SAMPLE_MAP"
echo " Dir   : $FASTQ_DIR"
echo "=================================================="

# ── main loop ────────────────────────────────────────────────────────────────
# awk handles both tab- and space-separated columns, and strips leading/trailing
# whitespace from each field, printing them tab-separated for safe read parsing.
while IFS=$'\t' read -r barcode sample; do

    # skip blank lines
    [[ -z "$barcode" || -z "$sample" ]] && continue

    for read_dir in R1 R2; do
        src="${FASTQ_DIR}/sample_${barcode}_${read_dir}.fastq.gz"
        dst="${FASTQ_DIR}/sample_${sample}_${read_dir}.fastq.gz"

        if [[ -f "$src" ]]; then
            # Warn if destination already exists (avoid silent overwrite)
            if [[ -f "$dst" ]]; then
                echo "SKIP (destination exists): $(basename "$dst")"
                ((skipped++))
                continue
            fi
            mv "$src" "$dst"
            echo "OK  : sample_${barcode}_${read_dir}.fastq.gz  ->  sample_${sample}_${read_dir}.fastq.gz"
            ((renamed++))
        else
            echo "WARN: not found -> sample_${barcode}_${read_dir}.fastq.gz"
            ((missing++))
        fi
    done

done < <(awk '{
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
    if ($1 != "" && $2 != "") printf "%s\t%s\n", $1, $2
}' "$SAMPLE_MAP")

# ── summary ──────────────────────────────────────────────────────────────────
echo "=================================================="
echo " Done."
echo "   Renamed  : $renamed files"
echo "   Missing  : $missing files (barcode not found in $FASTQ_DIR)"
echo "   Skipped  : $skipped files (destination already existed)"
echo "=================================================="

# Exit with error if any files were missing
(( missing > 0 )) && exit 1 || exit 0
