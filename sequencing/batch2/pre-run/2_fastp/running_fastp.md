## Checking the quality with fastp
fastp: Remove adapters, poly Q tails, merge reads (important for short frags). 

Fastp generates two different sets of files: unmerged and merged (the reads that could be merged, versus the ones that were long enough to be separated from each other). 


Parallelized version of the script.

```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=10G

# activate conda
module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/jsmontgomery/anaconda/fastp

WORKDIR=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed
cd $WORKDIR

run_fastp() {
    prefx=${1%_R1.fastq.gz}

    fastp \
        --in1 ${prefx}_R1.fastq.gz \
        --in2 ${prefx}_R2.fastq.gz \
        --out1 ${prefx}_R1.unmerged.fq.gz \
        --out2 ${prefx}_R2.unmerged.fq.gz \
        --merge \
        --merged_out ${prefx}.collapsed.fq.gz \
        --html ${prefx}.html \
        --json ${prefx}.json \
        --thread 2
}

export -f run_fastp

parallel -j $SLURM_NTASKS_PER_NODE run_fastp ::: *_R1.fastq.gz


```

## Analyzing fastp results

```
#!/bin/bash

module load jq

outfile="bacth2_fastp_summary.txt"

# Updated header (removed min/max)
echo -e "sample\tmean_len_before\tmean_len_after\tdup_rate\tinsert_peak\tq30_bases\tq30_percent\tgc_content\ttotal_reads" > "$outfile"

for json_file in *json; do

    sample=${json_file%.json}
    [ -f "$json_file" ] || continue

    jq -r --arg sample "$sample" '
        [
            $sample,
            .summary.before_filtering.read1_mean_length,
            .summary.after_filtering.read1_mean_length,
            .duplication.rate,
            (.insert_size.peak // "NA"),
            .summary.after_filtering.q30_bases,
            .summary.after_filtering.q30_rate,
            .summary.after_filtering.gc_content,
            .summary.after_filtering.total_reads
        ] | @tsv
    ' "$json_file" >> "$outfile"

done
```

### Calculate the number of base pairs per sample after fastp 

```
#!/bin/bash
#SBATCH --job-name=calc_bp
#SBATCH --output=calc_bp.out
#SBATCH --error=calc_bp.err
#SBATCH --time=10:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G   # memory per cpu-core

cd /scratch/midway2/rozennpineau/herbarium/batch2/demuxed

out=number_base_pairs_fastp.txt

# Create the output file only if it doesn't already exist
if [ ! -f "$out" ]; then
    echo -e "Sample\tnum_base_pairs" > "$out"
fi

for fq in *.fq.gz; do

    samp=$(basename "$fq")

    # Skip if this sample has already been processed
    if grep -q "^${samp}" "$out"; then
        echo "Skipping $samp (already processed)"
        continue
    fi

    echo "Processing $samp"

    num=$(zcat "$fq" | awk 'NR % 4 == 2 {sum += length($0)} END {print sum}')

    echo -e "$samp\t$num" >> "$out"

done
```