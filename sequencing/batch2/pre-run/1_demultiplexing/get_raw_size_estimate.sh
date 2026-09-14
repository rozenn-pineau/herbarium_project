#!/bin/bash
# sample_sizes.sh
# Sums R1 + R2 file sizes per sample, filters to a provided sample list,
# and writes a tab-separated summary (sample | size_GB) to an output file.
#
# Usage:
#   bash sample_sizes.sh <fastq_dir> <sample_list.txt> <output.txt> [col]
#
# Arguments:
#   fastq_dir       : directory containing sample_*_R1/R2.fastq.gz files
#   sample_list.txt : text file with sample IDs to keep (one per line,
#                     or a multi-column file — specify the column with [col])
#   output.txt      : path for the output summary file
#   col             : column number of sample IDs in sample_list.txt (default: 1)
#
# Example:
#   bash sample_sizes.sh /path/to/demux/ my_samples.txt sizes_summary.txt
#   bash sample_sizes.sh /path/to/demux/ metadata.txt   sizes_summary.txt 3

FASTQ_DIR="$1"
SAMPLE_LIST="$2"
OUTPUT="$3"
COL="${4:-1}"          # which column holds sample IDs (default: 1)

# ── input checks ─────────────────────────────────────────────────────────────
if [ -z "$FASTQ_DIR" ] || [ -z "$SAMPLE_LIST" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <fastq_dir> <sample_list.txt> <output.txt> [col]"
    exit 1
fi
if [ ! -d "$FASTQ_DIR" ]; then
    echo "ERROR: directory not found: $FASTQ_DIR"; exit 1
fi
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: sample list not found: $SAMPLE_LIST"; exit 1
fi

# ── step 1: sum R1 + R2 bytes per sample ─────────────────────────────────────
# ls -l gives exact byte counts (no human-rounding); awk splits filename on '_'
# to extract the sample number and accumulates bytes across R1 and R2.
ALL_SIZES=$(ls -l "${FASTQ_DIR}"/sample_*_R*.fastq.gz 2>/dev/null \
  | awk '
      NF >= 9 {
          filename = $NF
          # extract basename in case full path is printed
          n = split(filename, path, "/")
          base = path[n]
          # split on "_": sample _ <ID> _ R1/R2 .fastq.gz
          m = split(base, parts, "_")
          sample = parts[2]
          if (sample ~ /^[0-9]+$/) bytes[sample] += $5
      }
      END { for (s in bytes) print s, bytes[s] }
  ' | sort -n)

if [ -z "$ALL_SIZES" ]; then
    echo "ERROR: no sample_*_R*.fastq.gz files found in $FASTQ_DIR"
    exit 1
fi

echo "$ALL_SIZES" > /tmp/_sample_sizes_tmp.txt

# ── step 2: filter by sample list and convert bytes → GB ─────────────────────
found=0
missing=0

printf "sample\tsize_GB\n" > "$OUTPUT"

# Extract sample IDs from the requested column, skip blank lines and comments
while IFS= read -r line; do
    # skip blank lines and comment lines
    echo "$line" | grep -qE '^\s*(#|$)' && continue

    s=$(echo "$line" | awk -v col="$COL" '{print $col}' | tr -d '[:space:]')
    [ -z "$s" ] && continue

    match=$(grep "^${s} " /tmp/_sample_sizes_tmp.txt)
    if [ -n "$match" ]; then
        bytes=$(echo "$match" | awk '{print $2}')
        gb=$(awk "BEGIN {printf \"%.2f\", $bytes / 1073741824}")
        printf "%s\t%s\n" "$s" "$gb" >> "$OUTPUT"
        found=$((found + 1))
    else
        printf "%s\tNOT FOUND\n" "$s" >> "$OUTPUT"
        missing=$((missing + 1))
    fi
done < "$SAMPLE_LIST"

rm -f /tmp/_sample_sizes_tmp.txt

# ── summary ──────────────────────────────────────────────────────────────────
echo "=================================================="
echo " Output written to : $OUTPUT"
echo " Samples found     : $found"
echo " Samples NOT found : $missing  (listed as NOT FOUND in output)"
echo "=================================================="
