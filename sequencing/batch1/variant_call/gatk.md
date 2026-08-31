## Using Gatk to call variants 

We could use FreeBayes or gatk, but gatk can handle a larger number of samples.

Advantages of gatk: maximizing accuracy around indels or difficult genomic regions

It is a little more steps: 

- Run HaplotypeCaller independently on each BAM to produce a gVCF.
- Combine the gVCFs.
- Joint genotype all samples.


### Step 1. Prepare the reference (.fa + .fai + .dict)

```
module load gatk 

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta

#create sequence dictionary and fai index 
if [ ! -f ${ref%.fasta}.dict ]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R ${ref}
fi

if [ ! -f ${ref}.fai ]; then
    echo "Creating fasta index..."
    samtools faidx ${ref}
fi
```
This was very fast, the files were created in the reference directory directly. 

### Step 2. Split the reference genome into 100kb intervals. 

```
# ─── step 2: split genome into 100kb intervals 

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
workdir=/scratch/midway2/rozennpineau/herbarium/
intervaldir=/scratch/midway3/rozennpineau/herbarium/ref_intervals

# calculate number of 100kb windows from the fai index
n_windows=$(awk '{sum += $2} END {printf "%d", sum/100000}' ${ref}.fai)


echo "Splitting genome into 100kb windows..."
gatk SplitIntervals \
    -R ${ref} \
    -O ${intervaldir} \
    --scatter-count ${n_windows} \
    --subdivision-mode INTERVAL_SUBDIVISION \
    --interval-padding 0 


#parameters options explained: 
#--subdivision-mode INTERVAL_SUBDIVISION: simple scatter approach in which all output intervals have size equal to the total base count of the source list divided by the scatter count (except, possibly, in the last interval list).
#--interval-padding 0 : Amount of padding (in bp) to add to each interval you are including
#-L : genomic intervals over which to operate
#scatter count: number of output interval files to split into


# Count how many interval files were created
n_intervals=$(ls ${intervaldir}/*.interval_list | wc -l)
echo "[$(date)] Created ${n_intervals} interval files (~100kb each)"

```
This created 95 intervals (.list files). 

### Step 3. HaplotypeCaller - make GVCF files

```
# ── 1. define variables
bam_dir=/scratch/midway3/rozennpineau/herbarium/bams
work_dir=/scratch/midway3/rozennpineau/herbarium/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals
log_dir=${work_dir}/logs

mkdir -p ${work_dir}/tmp ${log_dir}

# ── 2. build list of all sample×interval combinations
combo_list=()
for bam in ${bam_dir}/*.dedup.sorted.bam; do
    for interval in ${interval_dir}/*.interval_list; do
        combo_list+=("${bam}:::${interval}")  # combine with a separator
    done
done

# ── 3. export variables
export work_dir ref log_dir

# ── 4. define the per-sample per-interval function
run_interval() {
    combo=$1
    bam="${combo%%:::*}"          # extract everything before :::
    interval_file="${combo##*:::}" # extract everything after :::

    sample=$(basename ${bam} .scaffolds.dedup.sorted.bam)
    interval_name=$(basename ${interval_file} .interval_list)
    out_dir=${work_dir}/${sample}
    mkdir -p ${out_dir}
    out_gvcf=${out_dir}/${sample}.${interval_name}.g.vcf.gz
    log_file=${log_dir}/${sample}.${interval_name}.log

    # skip if already done
    if [ -f ${out_gvcf} ]; then
        echo "[$(date)] Skipping ${sample} ${interval_name} — already exists"
        return 0
    fi

    echo "[$(date)] Running ${sample} on interval ${interval_name}..."
    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${bam} \
        -O ${out_gvcf} \
        -L ${interval_file} \
        -ERC BP_RESOLUTION \
        --tmp-dir ${work_dir}/tmp \
        --native-pair-hmm-threads 2 \
        --max-alternate-alleles 4 \
        > ${log_file} 2>&1


    if [ $? -eq 0 ]; then
        echo "[$(date)] Finished ${sample} ${interval_name}"
    else
        echo "[$(date)] ERROR on ${sample} ${interval_name}" >&2
        return 1
    fi
}

export -f run_interval

# ── 5. run parallel over all sample×interval combinations
printf '%s\n' "${combo_list[@]}" | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog ${log_dir}/parallel_joblog.txt \
             --halt soon,fail=1 \
             run_interval {}

echo "[$(date)] All jobs complete."


```


### Troubleshooting:


```
bam_dir=/scratch/midway3/rozennpineau/herbarium/bams
work_dir=/scratch/midway3/rozennpineau/herbarium/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals
log_dir=${work_dir}/logs
sample=$(basename ${bam} .scaffolds.dedup.sorted.bam)
interval_name=$(basename ${interval_file} .interval_list)
out_dir=${work_dir}/${sample}
mkdir -p ${out_dir}
out_gvcf=${out_dir}/${sample}.${interval_name}.g.vcf.gz


    gatk HaplotypeCaller \
        -R ${ref} \
        -I ${bam} \
        -O ${out_gvcf} \
        -L ${interval_file} \
        -ERC BP_RESOLUTION \
        --tmp-dir ${work_dir}/tmp \
        --native-pair-hmm-threads 2 \
        --max-alternate-alleles 4

# that worked !
```

Testing the parallel commands:
```
printf '%s\n' "${combo_list[@]}" | head -n 1 | \
    parallel --jobs 1 \
             --joblog ${log_dir}/parallel_joblog.txt \
             --halt soon,fail=1 \
             run_interval {}

#that worked !!

printf '%s\n' "${combo_list[@]}" | head -n 1 | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog /scratch/midway3/rozennpineau/herbarium/gvcf/logs/parallel_joblog.txt \
             --halt soon,fail=1 \
             run_interval {}

```


### Filters from Haplotype Caller that we could use to speed up process
--min-assembly-region-size 50 - minimum size of assembly region, default is 50, increase to skip very small regions

--max-assembly-region-size 300  maximum size of assembly region, default is 300, decrease to limit region complexity

--assembly-region-padding 100 - Number of additional bases of context to include around each assembly region, default 100, reduce to shrink flanking regions

--active-probability-threshold 0.002  # default 0.002, increase to call fewer active regions# how is this probability estimated in the first place?

--min-pruning 2 - Minimum support to not prune paths in the graph        default 2, increase to prune more aggressively - use with caution ! Using a prune factor of 1 (or below) will prevent any pruning from the graph, which is generally not ideal; it can make the calling much slower and even less accurate (because it can prevent effective merging of "tails" in the graph). Higher values tend to make the calling much faster, but also lowers the sensitivity of the results (because it ultimately requires higher depth to produce calls).


--max-unpruned-variants 100 - Maximum number of variants in graph the adaptive pruner will allow - # limit variants considered per active region, default is 100

--min-base-quality-score 10 - Minimum base quality required to consider a base for calling

--standard-min-confidence-threshold-for-calling 30 , minimum phred-scaled confidence threshold at which variants should be called, default is 30 (variant sites with QUAL equal or greater than)

--base-quality-score-threshold 18, base qualities below this threshold will be reduced to the minimum (6).

--max-alternate-alleles 6, Maximum number of alternate alleles to genotype, decreasing this number may speed up process

1. Prepare reference (.fa + .fai + .dict)       ← prerequisite for everything
2. Run SplitIntervals → interval shards          ← do this once, upfront
3. HaplotypeCaller (per sample × per shard)      ← parallelized
4. GenomicsDBImport (per shard)                  ← parallelized
5. GenotypeGVCFs (per shard)                     ← parallelized
6. GatherVcfs                                    ← merge back


gVCF = genomic VC --> records information about every position in the genome, not just the positions where a variant was found

Create a list of bam files (34 samples with more than 5x coverage):
```
ls /scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup/final_bams/for_variant_call/*bam > /scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup/final_bams/for_variant_call/bam_for_variant_call.list

```
Index and sort bam files (ran in an interative session):
```
for bam in *.scaffolds.dedup.bam; do
        prefx=${bam%.scaffolds.dedup.bam}
        sambamba sort -m 8GB --tmpdir tmp -t 2 -o ${prefx}.scaffolds.dedup.sorted.bam ${prefx}.scaffolds.dedup.bam
        sambamba index -t 2 ${prefx}.scaffolds.dedup.sorted.bam 
done
```

Make gVCF files:

```
#!/bin/bash
#SBATCH --job-name=gvcf
#SBATCH --output=gvcf.out
#SBATCH --error=gvcf.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=24G   # memory per cpu-core
#SBATCH --array=1-35

module load gatk
module load samtools

workdir=/scratch/midway2/rozennpineau/herbarium/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
bamlist=/scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup/final_bams/for_variant_call/bam_for_variant_call.list

cd $workdir

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bamlist)

sample=$(basename ${bam} .scaffolds.dedup.bam)

gatk --java-options "-Xmx20g" HaplotypeCaller \
    -R ${ref} \
    -I ${bam} \
    -O gvcf/${sample}.g.vcf.gz \
    -ERC GVCF \
    --native-pair-hmm-threads 8 #how many CPU threads to dedicate to the Pair-HMM calculations.

#The Pair-HMM calculates the probability that this read originated from each haplotype, taking into account: sequencing errors (base quality scores), insertions, deletions, mismatches
```