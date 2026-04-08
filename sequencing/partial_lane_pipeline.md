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


tool: fastp
use: remove adapters, poly Q tails, merge reads (important for short frags)
command line: herbarium_pipeline.sh sample ref_genome 
#sample is in the first positional argument in the bash script

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

### Map to reference genome with bwamem
```

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
# 2. Map merged (collapsed) reads to Reference Genome and calculate Endogenous DNA
bwa mem -t $threads -R "@RG\tID:$prefx\tSM:$prefx" $ref ${pathtobams}/${prefx}_1.fq.gz ${pathtobams}/${prefx}_2.fq.gz | samtools view -@ $threads -Sbh - >  /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/${prefx}.uns.bam
```

tool : bwamem
use: mapping
command line: 

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



