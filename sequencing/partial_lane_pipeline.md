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
#sambamba 1.0.1
#DeDup v0.12.9
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

### Calculate expected maximum coverage based on raw reads
Calculate the number of base pairs in unmerged and collapsed reads :
```
out=number_base_pairs_unmerged.txt
echo -e "Sample\tnum_base_pairs" > $out

for fq in */*.unmerged.fq.gz; do

  samp=$(basename "$fq")
  echo $samp
  num=$(zcat $fq | awk 'NR % 4 == 2 {sum += length($0)} END {print sum}')
  echo -e "$samp\t$num" >> $out

done
```

### Map to reference genome with bwamem - unmerged reads
```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
threads=6
cd /scratch/midway3/rozennpineau/herbarium_partial_lane/raw

for dir in ./* ; do
    cd "$dir" || continue

    for r1 in *_1.unmerged.fq.gz; do
        prefx=${r1%_1.unmerged.fq.gz}
        # 2. Map unmerged reads to Reference Genome, turn SAM to BAM 
        bwa mem -t $threads -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" $ref ${prefx}_1.unmerged.fq.gz ${prefx}_2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  /scratch/midway3/rozennpineau/herbarium_partial_lane/bams/${prefx}.unmerged.uns.bam

    done

    cd ..

done


#explanation: align (produces sam files) bwa mem -t (number of threads) -R (read group header line) $reference_genome reads_lane1 reads_lane2
#pipe: converts sam to bams with samtool view @ #threads -Sbh sam file input, convert to bam, with header 

```

The fastq header corresponds to : @Instrument:RunID:FlowcellID:Lane:Tile:X:Y Read:Filter:Control:Index

### Map to reference genome with bwamem - collapsed reads

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
threads=6
cd /scratch/midway3/rozennpineau/herbarium_partial_lane/raw

for dir in ./* ; do
    cd "$dir" || continue

    for r1 in *.collapsed.gz; do
        prefx=${r1%.collapsed.gz}
        # 2. Map merged (collapsed) reads to Reference Genome, turn SAM to BAM 
        bwa mem -t $threads -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" $ref ${prefx}.collapsed.gz | samtools view -@ $threads -Sbh - > /scratch/midway3/rozennpineau/herbarium_partial_lane/collapsed_bams/${prefx}.uns.bam

    done

    cd ..

done
```

### Count the number of base pairs in bams
Calculate the numebr of base pairs that were effectively mapped in both unmerged and collapased reads.

```
module load samtools
path=/scratch/midway3/rozennpineau/herbarium_partial_lane/merged_bams
cd $path
out=number_base_pairs_merged_bams.txt
echo -e "sample\tnum_base_pairs" > $out

for bam in *.renamed.final.sorted.bam; do

  name=${bam%.renamed.final.sorted.bam}
  echo "Calculating base pairs mapped in $name..."
  num=$(samtools stats $bam | grep "bases mapped (cigar):" | cut -f 3)
  echo -e "$name\t$num" >> $out

done


```

### Extract number of mapped reads from bams

```
outfile="number_mapped_reads_merged.txt"

echo -e  "sample\ttotal" > $outfile

threads=2
for file in *.bam; do
  name=$(basename "$file")
  echo "Counting mapped reads in $name"
  total=$(samtools view -@ $threads -c $file)
  echo -e "$name\t$total" >> $outfile
done

```

### Extract number of *un*mapped reads from bams

```
outfile="number_unmapped_reads.txt"

echo -e  "sample\ttotal" > $outfile

threads=2
for file in *.bam; do
  name=$(basename "$file")
  echo "Counting unmapped reads in $name"
  total=$(samtools view -@ $threads -c -f 4 $file)
  echo -e "$name\t$total" >> $outfile
done

```


### Sort bams
[Sambamba](https://lomereiter.github.io/sambamba/) is a high-performance, parallelized software tool written in the D programming language for fast processing of NGS SAM/BAM/CRAM alignment files.
```

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

threads=6
path_to_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/bams
cd $path_to_bams

for r1 in *.unmerged.uns.bam; do
        prefx=${r1%.unmerged.uns.bam}
        sambamba sort -m 15GB --tmpdir $path_to_bams/tmp -t $threads -o $path_to_bams/${prefx}.unmerged.sorted.bam $path_to_bams/${prefx}.unmerged.uns.bam
done

```



### Merge unmerged and collapsed bams

```
module load samtools
unmerged_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/bams/renamed/
collapsed_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/collapsed_bams
output_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/merged_bams

cd $unmerged_bams
for r1 in *.unmerged.sorted.bam; do

        prefx=${r1%.unmerged.sorted.bam}
        samtools merge $output_folder/${prefx}.final.sorted.bam $unmerged_bams/${prefx}.unmerged.sorted.bam $collapsed_bams/${prefx}.sorted.bam

done

```

Sample herb189 was run on two different lanes, merge the bams now:

```
#merge

samtools merge -f herb189.final.sorted.bam herb189_CKDL260004898-1A_23FFGFLT4_L2.final.sorted.bam herb189_CKDL260004898-1A_23F5GKLT4_L1.final.sorted.bam

````

Command lines for the files that did not work in the loop (herb5). 

```

#check the header
samtools view herb5_CKDL260004894-1A_23F5GKLT4_L1_1.sorted.bam | head

#merge
unmerged_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/bams
collapsed_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/collapsed_bams
output_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams

cd $unmerged_bams
samtools merge -f $output_folder/herb5_CKDL260004894-1A_23F5GKLT4_L1.final.sorted.bam $unmerged_bams/herb5_CKDL260004894-1A_23F5GKLT4_L1.unmerged.sorted.bam $collapsed_bams/herb5_CKDL260004894-1A_23F5GKLT4_L1.sorted.bam


````


### Calculate mapped reads after merging for comparison
Using samtools_num_reads.sh script.

### Calculate endogenous DNA
echo "EndogenousDNA" >> ${prefx}.log
python -c "print(float($mapped)/ $total)" >> ${prefx}.log

### Map Damage
[MapDamage](https://ginolhac.github.io/mapDamage/): calculate deamination levels per sample, across the reads, and rescale per-base quality.

I had some problems with dependencies when using my normal conda environment, so I created another one, called "mapdam".
```
#!/bin/bash
#SBATCH --job-name=mapdamage
#SBATCH --output=mmapdamage.out
#SBATCH --error=mapdamage.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate mapdam


ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
cd /scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams

for r1 in *.final.sorted.bam; do

        prefx=${r1%.final.sorted.bam}
        mapDamage -i ${prefx}.final.sorted.bam -r $ref --rescale

done

#The --rescale parameter can be optionally used to rescale quality scores of likely damaged positions in the reads. A new BAM file is constructed by downscaling quality values for misincorporations likely due to ancient DNA damage according to their initial qualities, position in reads and damage patterns.
```

### Estimate duplication rates

Trying with two different tools, DeDup and Picard. 

With DeDup:

I first have to split the bam into chromosomes, as it keeps running out of memory on the full bam file. 

```
#create a subfolder for each file
#split bam into chromosomes

!!! There is a problem with this script in that it looks at more and more bams as it goes, copied the "herb" files into the scaffold folders.

Make sure to test it before running again.

it might be in the chr loop


bam=herb5_CKDL260004894-1A_23F5GKLT4_L1.renamed.final.sorted.bam

name=${bam%.renamed.final.sorted.bam}
mkdir $name
for chr in $(samtools idxstats $bam | cut -f1); do
    echo "Splitting $chr..."
    samtools view -b $bam "$chr" > $name/"${chr}.bam" #split into scaffolds
    samtools index $name/"${chr}.bam" #index
done

```

[DeDup](https://github.com/apeltzer/DeDup/blob/master/README.md) tool is a PCR duplicate removal tool of paired-end and merged sequenced data designed for ultra-short DNA (e.g. ancient DNA).

DeDup looks for reads with the same start and end coordinates, and whether they have exactly the same sequence. The main difference of DeDup versus e.g. samtools markduplicates is that DeDup considers both ends of a read, not just the start position, so it is more precise in removing actual duplicates without penalising often already low aDNA data. (from: https://nf-co.re/eager/2.4.1/docs/output#dedup)


Run DeDup:
```
module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
conda activate mapdam

bam_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams
cd $bam_folder
output_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams/dedup_marked_dups

for bam in *.renamed.final.sorted.bam; do
        name=${bam%.renamed.final.sorted.bam}
        mkdir -p $output_folder/$name
        dedup -i $bam -o $output_folder/$name
done
```

Initially, the run was not able to write the files to the folder because of a memory error. I increased the memory allocation and started the run again. 

With picard: 
```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

bam_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams
cd $bam_folder
output_folder=/scratch/midway3/rozennpineau/herbarium_partial_lane/final_bams/picard_marked_dups

for bam in .renamed.final.sorted.bam; do
        name=${bam%.renamed.final.sorted.bam}
        picard MarkDuplicates I=${name}.renamed.final.sorted.bam O=$output_folder/${name}.dup.renamed.final.sorted.bam M=$output_folder/${name}.dup.txt
done

```

### Retrieve duplication rate from files after running DeDup

```
echo -e "sample\tmean duplication rate" > mean_dup_rate.txt
for sample_dir in herb*/; do
cd $sample_dir
  echo -e "scaffold\tduplication rate" > dup_rate_summary.txt
  for scaffold_dir in [Ss]*/; do
    scaffold_name=$(basename "$scaffold_dir")
    
    dup_rate=$(awk '/Duplication Rate/ {print $3}' "$scaffold_dir/$scaffold_name.log") #Retrieve duplication rate

    echo -e "${scaffold_name}\t${dup_rate}" >> dup_rate_summary.txt
  done
  #calculate mean over all scaffolds
  mean=$(awk 'NR>1 {sum += $2; n++} END {print sum/n}' dup_rate_summary.txt)
cd ..
echo -e "$sample_dir\t$mean" >> mean_dup_rate.txt
done
```

Calculate mean duplication rate for Scaffolds 1-16:
```
echo -e "Sample\tScaffold1-16_duplication_rate" > scaffold1-16_mean_dup_rate.txt

for dir in ./herb*; do
  #echo $dir
  dup_rate=$(cat $dir/dup_rate_summary.txt | grep Scaffold | awk 'NR>1 {sum += $2; n++} END {print sum/n}')
  echo -e "${dir}\t${dup_rate}" >> scaffold1-16_mean_dup_rate.txt
done
```

### Merge files after running DeDup
Dedup create a subfolder for each file, and I had to scaffold the samples to run it, so I am left with 1000 + folders for each file. I would like to merge the deduped bams to calculate samtools depth. 

```
start_dir=/scratch/midway3/rozennpineau/herbarium_partial_lane/merged_bams/scaffolds

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bamtools

ulimit -n 4096 # increase upper limit of number of files that can be opened at once

cd $start_dir
for dir in ./*; do
    cd $dir

    realpath */[s]*rmdup.bam > bams_to_merge.list
    bamtools merge -list bams_to_merge.list -out $dir.scaffolds.dedup.bam

    cd ..
done
```


### Calculate depth at each step
With samtools:
```
module load samtools

out="average_bam_depths.txt"
echo -e "Sample\tAverage_Depth" > $out

for bam in *.dup.renamed.final.sorted.bam; do
    name=${bam%.dup.renamed.final.sorted.bam}
    AVG_DEPTH=$(samtools depth "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    echo -e "${name}\t${AVG_DEPTH}" >> $out
done
```





