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

module load samtools

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
cd /scratch/midway3/rozennpineau/herbarium_partial_lane/raw

for dir in ./* ; do
    cd "$dir" || continue

    for r1 in *_1.unmerged.fq.gz; do
        prefx=${r1%_1.unmerged.fq.gz}
        # 2. Map merged (collapsed) reads to Reference Genome, turn SAM to BAM 
        bwa mem -t 6 -R "@RG\tID:$NA\tSM:$NA\tPL:ILLUMINA\tLB:$NA" $ref ${prefx}_1.unmerged.fq.gz ${prefx}_2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  /scratch/midway3/rozennpineau/herbarium_partial_lane/bams/${prefx}.unmerged.uns.bam

    done

    cd ..

done



#explanation: align (produces sam files) bwa mem -t (number of threads) -R (read group header line) $reference_genome reads_lane1 reads_lane2
#pipe: converts sam to bams with samtool view @ #threads -Sbh sam file input, convert to bam, with header 

```
Troubleshooting code:

```
for r1 in *_L1.short.unmerged.fq.gz ; do
        prefx=${r1%_L1.short.unmerged.fq.gz}
        # 2. Map merged (collapsed) reads to Reference Genome, turn SAM to BAM 
        bwa mem -t 6 -R "@RG\tID:$prefx\tSM:$prefx\tPL:ILLUMINA\tLB:$prefx" $ref ${prefx}_L1.short.unmerged.fq.gz ${prefx}_L2.short.unmerged.fq.gz | samtools view -@ 6 -Sbh - > ${prefx}.unmerged.uns.bam

done
```
The fastq header corresponds to : @Instrument:RunID:FlowcellID:Lane:Tile:X:Y Read:Filter:Control:Index



### Extract results from fastp files

Using the *json* block file to extract the information for each fastp report:

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
### Extract some stats from bams

```
threads=2
for file in *.bam; do
  name=$(basename "$file")
  total=$(samtools view -@ $threads -c $file)
  echo -e "TotalReads\n$total" >> $file.log
done
```
sambamba sort -m 15GB --tmpdir $path/bams/tmp -t $threads -o /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/$prefx.sorted.scaled.bam /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/${prefx}.scaled.bam
samtools merge ${prefx}.final.sorted.bam ${prefx}.unmerg.sorted.bam ${prefx}.merged.sorted.bam
rm ${prefx}.uns.bam
mapped=$(samtools view -@ $threads -c $prefx.bam)
echo -e "MappedReads\n$mapped" >> ${prefx}.log
echo "EndogenousDNA" >> ${prefx}.log
python -c "print(float($mapped)/ $total)" >> ${prefx}.log



