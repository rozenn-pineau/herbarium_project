#Goal: to record the pipeline for the partial lane sequencing. Draw from Julia's GitHub + Science paper.

Steps should include:
Comparison of C to T deamination between Julia's raw reads and the new raw reads using MapDamage. 
  --> 	does the lib prep decrease deamination?

Comparison before and after merging reads
  --> map before and after merging reads: set of reads paired and merged, map twice, merge resulting bams, compare results

Comparison of the number of duplicates between sets of herbarium sequence data AND between sample quality/concentrations. 

Compare mean fragment size between between sets of herbarium sequence data AND between sample quality/concentrations. 

### Software versions
```
#fastp version 0.23.4
#bwamem2-2.2.1
```

### Check that all files were downloaded correctly from Novogene
An MD5 checksum is a 32-character hexadecimal number that is computed on
a file. If two files have the same MD5 checksum value, it is highly probable that
the two files are the same. It is generally used to check data integrity.

```
#check the md5 string, printed to the screen ("OK")
for dir in ./*
do
    cd $dir
    md5sum -c MD5.txt 
    cd ../
done
#all samples looked good.
```

### Check quality with fastp


Remove adapters, poly Q tails, merge reads (important for short frags).
command line: herbarium_pipeline.sh sample ref_genome 
#sample is in the first positional argument in the bash script

Note: the “unmerged” files do not include reads that were successfully merged. Fastp generates two different sets of files: unmerged and merged (the reads that could be merged, versus the ones that were long enough to be separated from each other). I will do two alignments, merged and unmerged and comapre outputs for both.

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner


cd /scratch/midway2/rozennpineau/herbarium_partial_lane/raw

for dir in ./* ; do
    cd "$dir" || continue

    for r1 in *_1.fq.gz; do
        prefx=${r1%_1.fq.gz}

        fastp \
            --in1 ${prefx}_1.fq.gz \
            --in2 ${prefx}_2.fq.gz \
            --out1 ${prefx}_1.unmerged.fq.gz \
            --out2 ${prefx}_2.unmerged.fq.gz \
            --merge \
            --merged_out ${prefx}.collapsed.gz
    done

    cd ..
done
```
### Analyze fastp results
Using the *json* block file to extract the information for each fastp report:

```
#!/bin/bash

outfile="fastp_summary.txt"

# Updated header (removed min/max)
echo -e "sample\tmean_len_before\tmean_len_after\tdup_rate\tinsert_peak\tq30_bases\tq30_percent\tgc_content\ttotal_reads" > "$outfile"

for dir in */; do
    sample=$(basename "$dir")
    json_file="${dir}/fastp.json"

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



### Map to reference genome with bwamem
```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
cd /scratch/midway2/rozennpineau/herbarium_partial_lane/raw

for dir in ./* ; do
    cd "$dir" || continue

    for r1 in *_1.unmerged.fq.gz; do
        prefx=${r1%_1.unmerged.fq.gz}
        # 2. Map merged (collapsed) reads to Reference Genome, turn SAM to BAM 
        bwa mem -t 6 -R "@RG\tID:$prefx\tSM:$prefx" $ref ${prefx}_1.unmerged.fq.gz ${prefx}_2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  /scratch/midway2/rozennpineau/herbarium_partial_lane/bams/${prefx}.unmerged.uns.bam

    done

    cd ..

done

```

Custom scripts for the two remaining files that did ot make it within the 10 hours alloted for the job:

```
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
cd /scratch/midway2/rozennpineau/herbarium_partial_lane/raw

bwa mem -t 6 -R "@RG\tID:$prefx\tSM:$prefx" $ref herb5_CKDL260004894-1A_23F5GKLT4_L1_1.unmerged.fq.gz herb5_CKDL260004894-1A_23F5GKLT4_L1_2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  /scratch/midway2/rozennpineau/herbarium_partial_lane/bams/herb5_CKDL260004894-1A_23F5GKLT4_L1_1.unmerged.uns.bam

bwa mem -t 6 -R "@RG\tID:$prefx\tSM:$prefx" $ref herb449_CKDL260004902-1A_23FFGFLT4_L6_1.unmerged.fq.gz herb449_CKDL260004902-1A_23FFGFLT4_L6_2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  /scratch/midway2/rozennpineau/herbarium_partial_lane/bams/herb449_CKDL260004902-1A_23FFGFLT4_L6_1.unmerged.uns.bam

```

### fastp files

Using tje *json* block file to extract the information for each fastp report:

```

#!/bin/bash

# Output file
outfile="fastp_summary.txt"

# Header
echo -e "sample\tmean_len_before\tmean_len_after\tdup_rate\tinsert_peak\tinsert_min\tinsert_max\tq30_bases\tq30_percent\tgc_content\ttotal_reads" > $outfile

# Loop through directories
for dir in */; do
    sample=$(basename "$dir")
    html="${dir}/fastp.html"

    # Skip if no html file
    [ -f "$html" ] || continue

    # Extract JSON block from HTML
    json=$(sed -n '/<script id="json" type="application\/json">/,/<\/script>/p' "$html" | sed '1d;$d')

    # Parse values using jq
    echo "$json" | jq -r --arg sample "$sample" '
        [
            $sample,
            .summary.before_filtering.read1_mean_length,
            .summary.after_filtering.read1_mean_length,
            .duplication.rate,
            .insert_size.peak,
            .insert_size.min,
            .insert_size.max,
            .summary.after_filtering.q30_bases,
            .summary.after_filtering.q30_rate,
            .summary.after_filtering.gc_content,
            .summary.after_filtering.total_reads
        ] | @tsv
    ' >> $outfile

done

```


download to laptop to open in Chrome
```
scp -r rozennpineau@midway3.rcc.uchicago.edu:/scratch/midway3/rozennpineau/herbarium_partial_lane/raw/\*/\*.html /Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My\ Drive/Work/9.Science/4.Herbarium/4.Sequencing/4.PartialLane/fastp


```
tool : DeDup
use: deduplication
command line: 

tool : MapDamage https://ginolhac.github.io/mapDamage/
use: quantify damage patterns (C to T substitutions because of C deamination) --> use to rescale per base quality score
command line: 
(check relationship with sample age)

tool : FreeBayes (maybe GTAK for 600+ samples?)
use: deduplication
command line: 
--use-best-n-alleles 4, --report-monomorphic in 100kb regions



