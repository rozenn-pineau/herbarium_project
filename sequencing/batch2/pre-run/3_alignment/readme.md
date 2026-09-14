


Align reads to reference genome using baw mem.

Align unmerged and collapsed reads separately.

### Aligning unmerged reads

Script parallelized with GNU parallel: 
```
#!/bin/bash
#SBATCH --job-name=bwa_unmerged
#SBATCH --output=bwa_unmerged.out
#SBATCH --error=bwa_unmerged.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7G

module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
module load samtools
module load parallel

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged
threads=$SLURM_CPUS_PER_TASK
working_dir=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed

cd $working_dir

mkdir -p $out/tmp

BWA_ENV=/project/kreiner/rpineau/bwa
SAMBAMBA_ENV=/project/kreiner/rpineau/sambamba

process_sample() {
    r1=$1
    prefx=${r1%_R1.unmerged.fq.gz}
    output=$out/${prefx}.unmerged.uns.bam

    if [ -s "$output" ]; then
        echo "$output already exists, skipping."
    else
        echo "Processing $prefx..."

        # Map unmerged reads
        $BWA_ENV/bin/bwa mem \
            -t $threads \
            -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" \
            $ref \
            ${prefx}_R1.unmerged.fq.gz \
            ${prefx}_R2.unmerged.fq.gz \
            | samtools view -@ $threads -Sbh - > $output

        # Sort bams
        $SAMBAMBA_ENV/bin/sambamba sort -m 8GB --tmpdir $out/tmp -t $threads -o $out/${prefx}.unmerged.sorted.bam $output

    fi
}

export -f process_sample
export out ref threads BWA_ENV SAMBAMBA_ENV

cd $working_dir
parallel -j $SLURM_NTASKS_PER_NODE process_sample ::: *_R1.unmerged.fq.gz

# for testing: parallel -j 1 process_sample sample_102_R1.unmerged.fq.gz
```

### Find samples that need alignment - unmerged

```
# Check each sample and keep only those that need processing
cd /scratch/midway2/rozennpineau/herbarium/batch2/demuxed

for r1 in *_R1.unmerged.fq.gz; do
    prefx=${r1%_R1.unmerged.fq.gz}
    output=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged/${prefx}.unmerged.uns.bam
    
    # If quickcheck fails (file missing or corrupt), add to todo list
    if ! samtools quickcheck $output 2>/dev/null; then
        echo $r1
    fi
done > unmerged_to_do.txt
```
### Job array - unmerged
I was not using resources efficiently enough, so I am now trying job arrays instead of GNU-parallel:

```
#!/bin/bash
#SBATCH --job-name=bwa_unmerged
#SBATCH --output=logs/bwa_%a.out
#SBATCH --error=logs/bwa_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7G

module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
module load samtools

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged
threads=$SLURM_CPUS_PER_TASK
working_dir=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed

BWA_ENV=/project/kreiner/rpineau/bwa
SAMBAMBA_ENV=/project/kreiner/rpineau/sambamba

mkdir -p $out/tmp logs

# Pick the Nth sample based on array task ID
cd $working_dir
r1=$(ls *_R1.unmerged.fq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
prefx=${r1%_R1.unmerged.fq.gz}
output=$out/${prefx}.unmerged.uns.bam

if samtools quickcheck $output 2>/dev/null; then
    echo "$output already exists and is valid, skipping."
else
    echo "Processing $prefx..."

    $BWA_ENV/bin/bwa mem \
        -t $threads \
        -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" \
        $ref \
        ${prefx}_R1.unmerged.fq.gz \
        ${prefx}_R2.unmerged.fq.gz \
        | samtools view -@ $threads -Sbh - > $output

    $SAMBAMBA_ENV/bin/sambamba sort \
        -m 15GB \
        --tmpdir $out/tmp \
        -t $threads \
        -o $out/${prefx}.unmerged.sorted.bam \
        $output
fi
```
Submit the batch job with 

```
# Submit with the exact right array size
N=$(wc -l < /scratch/midway2/rozennpineau/herbarium/batch2/demuxed/unmerged_to_do.txt)
sbatch --array=1-$N%8 run_bwa_unmerged_array.sh #max 8 nodes at one time
```
### Find samples that need alignment - collapsed

```
# Check each sample and keep only those that need processing
cd /scratch/midway2/rozennpineau/herbarium/batch2/demuxed

for r1 in *.collapsed.fq.gz; do
    prefx=${r1%.collapsed.fq.gz}
    output=/scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed/${prefx}.collapsed.uns.bam
    
    # If quickcheck fails (file missing or corrupt), add to todo list
    if ! samtools quickcheck $output 2>/dev/null; then
        echo $r1
    fi
done > collapsed_to_do.txt

```
### Aligning collapsed reads

```
```
#!/bin/bash
#SBATCH --job-name=bwa_unmerged
#SBATCH --output=logs/bwa_%a.out
#SBATCH --error=logs/bwa_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7G

module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
module load samtools

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged
threads=$SLURM_CPUS_PER_TASK
working_dir=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed

BWA_ENV=/project/kreiner/rpineau/bwa
SAMBAMBA_ENV=/project/kreiner/rpineau/sambamba

mkdir -p $out/tmp logs

# Pick the Nth sample based on array task ID
cd $working_dir
r1=$(ls *_R1.unmerged.fq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
prefx=${r1%_R1.unmerged.fq.gz}
output=$out/${prefx}.unmerged.uns.bam

if samtools quickcheck $output 2>/dev/null; then
    echo "$output already exists and is valid, skipping."
else
    echo "Processing $prefx..."

    $BWA_ENV/bin/bwa mem \
        -t $threads \
        -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" \
        $ref \
        ${prefx}.collapsed.fq.gz \
        | samtools view -@ $threads -Sbh - > $output

    $SAMBAMBA_ENV/bin/sambamba sort \
        -m 15GB \
        --tmpdir $out/tmp \
        -t $threads \
        -o $out/${prefx}.unmerged.sorted.bam \
        $output
fi
```
Submit the batch job with 

```
# Submit with the exact right array size
N=$(wc -l < /scratch/midway2/rozennpineau/herbarium/batch2/demuxed/collapsed_to_do.txt)
sbatch --array=1-$N%8 run_bwa_unmerged_array.sh #max 8 nodes at one time
```


#SBATCH --job-name=bwa_collapsed
#SBATCH --output=bwa_collapsed.out
#SBATCH --error=bwa_collapsed.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7G

module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh

#module load python/anaconda-2022.05
#source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh

module load samtools
module load parallel

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed
threads=$SLURM_CPUS_PER_TASK
working_dir=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed

mkdir -p $out/tmp

BWA_ENV=/project/kreiner/rpineau/bwa
SAMBAMBA_ENV=/project/kreiner/rpineau/sambamba

process_sample() {
    r1=$1
    prefx=${r1%.collapsed.fq.gz}
    output=$out/${prefx}.collapsed.uns.bam

    if [ -s "$output" ]; then
        echo "$output already exists, skipping."
    else
        echo "Processing $prefx..."

        # Map unmerged reads
        $BWA_ENV/bin/bwa mem \
            -t $threads \
            -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" \
            $ref \
            ${prefx}.collapsed.fq.gz \
            | samtools view -@ $threads -Sbh - > $output

        # Sort bams
        $SAMBAMBA_ENV/bin/sambamba sort -m 7GB --tmpdir $out/tmp -t $threads -o $out/${prefx}.collapsed.sorted.bam $output

    fi
}

export -f process_sample
export out ref threads BWA_ENV SAMBAMBA_ENV

cd $working_dir
parallel -j $SLURM_NTASKS_PER_NODE process_sample ::: *.collapsed.fq.gz
```

### Re-align 
Some files that were corrupted because the run ended early.
(1) find them
(2) re-align and sort


## Notes from this run on 532 (+4) samples 
Some samples were aligned but were not sorted (sambamba) so I had to go back and sort them. 

It took 4 days for those sampels to align (initially with 4 nodes instead of 8, going much faster with 8/10 nodes). 

The next step is to de-deduplicate the reads (new readme). 