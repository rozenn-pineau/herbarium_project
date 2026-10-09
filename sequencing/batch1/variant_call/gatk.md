## Using Gatk to call variants 

We could use FreeBayes or gatk, but gatk can handle a larger number of samples.

Advantages of gatk: maximizing accuracy around indels or difficult genomic regions

It is a little more steps: 

- Run HaplotypeCaller independently on each BAM to produce a gVCF.
- Combine the gVCFs.
- Joint genotype all samples.


### Step 1. Prepare the reference (.fa + .fai + .dict)

See the readme in 2022_dataset for this step.

### Step 2. Split the reference genome into 5 MB intervals. 

See the readme in 2022_dataset for this step.

### Step 3. HaplotypeCaller - make GVCF files

```
#!/bin/bash
#SBATCH --job-name=make_gvcfs
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=7:00:00
#SBATCH --array=1-1%40     # one task per BAM, 40 running at once
#SBATCH --output=logs/make_gvcfs_%A_%a.out
#SBATCH --error=logs/make_gvcfs_%A_%a.err

module load gatk samtools parallel

#load directories
main_dir=/scratch/midway3/rozennpineau/herbarium/batch1
bam_dir=${main_dir}/bams
work_dir=${main_dir}/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
log_dir=${work_dir}/logs
mkdir -p ${work_dir}/tmp ${log_dir}

# one sample per array task
bam=$(ls ${bam_dir}/*.scaffolds.dedup.sorted.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(basename ${bam} .scaffolds.dedup.sorted.bam)
out_dir=${work_dir}/${sample}
mkdir -p ${out_dir}

run_interval() {
    interval_file=$1
    interval_name=$(basename ${interval_file} .interval_list)
    out_gvcf=${out_dir}/${sample}.${interval_name}.g.vcf.gz
    log_file=${log_dir}/${sample}.${interval_name}.log

    # done only if the index exists (written last)
    [ -f ${out_gvcf}.tbi ] && return 0
    rm -f ${out_gvcf} ${out_gvcf}.tbi

    gatk --java-options "-Xmx5g" HaplotypeCaller \
        -R ${ref} -I ${bam} -O ${out_gvcf} -L ${interval_file} \
        -ERC GVCF \
        --tmp-dir ${work_dir}/tmp \
        --native-pair-hmm-threads 1 \
        --max-alternate-alleles 4 \
        > ${log_file} 2>&1
}
export -f run_interval
export work_dir ref log_dir bam sample out_dir

parallel --jobs ${SLURM_CPUS_PER_TASK} --retries 2 \
         --joblog ${log_dir}/joblog_${sample}.txt \
         run_interval ::: ${interval_dir}/*.interval_list

#compare number of shards expected versus obtained and report error
n_ok=$(ls ${out_dir}/*.g.vcf.gz.tbi | wc -l)
n_expected=$(ls ${interval_dir}/*.interval_list | wc -l)
[ "$n_ok" -eq "$n_expected" ] || { echo "ERROR: ${n_ok}/${n_expected} shards done" >&2; exit 1; }

echo "[$(date)] ${sample} complete."
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
### Time it took
Job #1 on 34 samples split into 6000+ 100kb windows separated into 94 intervals, 9 jobs at a time: 
15:46:04 to 20:29:10 --> 4 hours, 43 minutes with option BP_RESOLUTION
Job #2
11:28 to 16:04 --> 4 hours and 36 minutes with GVCF option

Difference in file sizes
2.2G	unmerged_gvcfs/ #BP_RESOLUTION
885M	unmerged_gvcfs/ #BP_RESOLUTION


### Step 4 - Merging gvcf files
Maybe add this step to the previous script !
```
#!/bin/bash
#SBATCH --job-name=gatk_merge
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --output=gatk_merge.out
#SBATCH --error=gatk_merge.err

# modules
module load gatk
module load parallel

# ── 1. define variables
work_dir=/scratch/midway3/rozennpineau/herbarium/gvcf_BP_RESOLUTION
unmerged_dir=${work_dir}/unmerged_gvcfs
merged_dir=${work_dir}/merged_gvcfs
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals
log_dir=${merged_dir}/logs
db_dir=${work_dir}/genomicsdb   # output directory for GenomicsDBImport

mkdir -p ${log_dir} ${db_dir} ${merged_dir}

echo "Running with ${SLURM_CPUS_PER_TASK} parallel jobs"

# ──────────────────────────────────────────────────────────────────────────────
# STEP 2: GatherVcfs — merge per-interval GVCFs into one GVCF per sample
# ──────────────────────────────────────────────────────────────────────────────

echo "[$(date)] Starting Step 2: GatherVcfs per sample..."

gather_sample() {
    sample_dir=$1
    sample=$(basename ${sample_dir})
    out_gvcf=${merged_dir}/${sample}.merged.g.vcf.gz
    log_file=${log_dir}/${sample}.gather.log

    # skip if already done
    if [ -f ${out_gvcf} ]; then
        echo "[$(date)] Skipping ${sample} — merged GVCF already exists"
        return 0
    fi

    # build sorted list of input GVCFs for this sample
    # sorting numerically ensures intervals are in the correct genome order
    input_vcfs=$(ls ${sample_dir}/${sample}.*.g.vcf.gz | sort -V | \
                 awk '{print "--INPUT "$1}' | tr '\n' ' ')

    echo "[$(date)] Gathering ${sample}..."
    gatk GatherVcfs \
        ${input_vcfs} \
        --OUTPUT ${out_gvcf} \
        > ${log_file} 2>&1

    # index the merged GVCF
    gatk IndexFeatureFile \
        -I ${out_gvcf} \
        >> ${log_file} 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date)] Finished gathering ${sample}"
    else
        echo "[$(date)] ERROR gathering ${sample}" >&2
        return 1
    fi
}

export -f gather_sample
export work_dir log_dir merged_dir unmerged_dir

# run GatherVcfs in parallel across samples
ls -d ${unmerged_dir}/sample*/ | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog ${log_dir}/gather_joblog.txt \
             --halt soon,fail=1 \
             gather_sample {}

```

This step worked and was very fast - add it at the end of the previous script?

### Step 5 - GenomicsDBImport
What does [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) do: 
Joint genotyping with GenotypeGVCFs needs to consider all samples simultaneously, but reading a lot (500) GVCFs directly is slow and inefficient. GenomicsDBImport solves this by consolidating them into an optimized database format. It reorganizes the data from being sample-oriented to site-oriented, to then be able to call variantas across all samples at each position in the genome.


Options for GenomicsDBImport that we night want/need to play with:
--batch-size 0	Batch size controls the number of samples for which readers are open at once and therefore provides a way to minimize memory consumption. However, it can take longer to complete. Use the consolidate flag if more than a hundred batches were used. This will improve feature read time. batchSize=0 means no batching (i.e. readers for all samples will be opened at once) Defaults to 0.

--reader-threads
1	How many simultaneous threads to use when opening VCFs in batches; higher values may improve performance when network latency is an issue

```
#!/bin/bash
#SBATCH --job-name=gatk_gdbimport
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --time=36:00:00
#SBATCH --output=gatk_gdbimport.out
#SBATCH --error=gatk_gdbimport.err

# modules
module load gatk
module load parallel

# ── 1. define variables
work_dir=/scratch/midway3/rozennpineau/herbarium/gvcf_BP_RESOLUTION
merged_dir=${work_dir}/merged_gvcfs
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals
log_dir=${work_dir}/logs
db_dir=${work_dir}/genomicsdb   # output directory for GenomicsDBImport

mkdir -p ${work_dir}/tmp

# ── 2. build the sample map: two columns — sample name and path to merged GVCF
# GenomicsDBImport requires this file
sample_map=${work_dir}/sample_map.txt
> ${sample_map}   # empty the file if it exists

for gvcf in ${merged_dir}/*.merged.g.vcf.gz; do
    sample=$(basename ${gvcf} .merged.g.vcf.gz)
    echo -e "${sample}\t${gvcf}" >> ${sample_map}
done

echo "[$(date)] Sample map written to ${sample_map} with $(wc -l < ${sample_map}) samples"

# ── 3. define per-interval GenomicsDBImport function
run_genomicsdb() {
    interval_file=$1
    interval_name=$(basename ${interval_file} .interval_list)
    db_path=${db_dir}/${interval_name}
    log_file=${log_dir}/${interval_name}.genomicsdb.log

    # skip if already done
    if [ -d ${db_path} ]; then
        echo "[$(date)] Skipping ${interval_name} — database already exists"
        return 0
    fi

    echo "[$(date)] Running GenomicsDBImport on interval ${interval_name}..."
    gatk GenomicsDBImport \
        --sample-name-map ${sample_map} \
        --genomicsdb-workspace-path ${db_path} \
        -L ${interval_file} \
        --reader-threads 2 \
        --batch-size 50 \
        --tmp-dir ${work_dir}/tmp \
        > ${log_file} 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date)] Finished GenomicsDBImport for interval ${interval_name}"
    else
        echo "[$(date)] ERROR on GenomicsDBImport for interval ${interval_name}" >&2
        return 1
    fi
}

export -f run_genomicsdb
export db_dir sample_map work_dir log_dir

# ── 4. run GenomicsDBImport in parallel across intervals
ls ${interval_dir}/*.interval_list | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog ${log_dir}/genomicsdb_joblog.txt \
             --halt soon,fail=1 \
             run_genomicsdb {}


```
### Step 6 - Genotype GVCFs

[GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/21905118377755-GenotypeGVCFs)

```
#!/bin/bash
#SBATCH --job-name=genotypegvcfs
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16G
#SBATCH --time=36:00:00
#SBATCH --output=gatk_genotypegvcfs.out
#SBATCH --error=gatk_genotypegvcfs.err

# modules
module load gatk
module load parallel

# ── 1. define variables
work_dir=/scratch/midway3/rozennpineau/herbarium/gvcf_BP_RESOLUTION
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals
db_dir=${work_dir}/genomicsdb        # GenomicsDBImport output from Step 3
genotype_dir=${work_dir}/genotyped   # output directory for this step
log_dir=${work_dir}/logs

mkdir -p ${genotype_dir} ${log_dir} ${work_dir}/tmp

echo "Running with ${SLURM_CPUS_PER_TASK} parallel jobs"

# ── 2. define per-interval GenotypeGVCFs function
run_genotype() {
    interval_file=$1
    interval_name=$(basename ${interval_file} .interval_list)
    db_path=${db_dir}/${interval_name}
    out_vcf=${genotype_dir}/${interval_name}.genotyped.vcf.gz
    log_file=${log_dir}/${interval_name}.genotype.log

    # skip if already done
    if [ -f ${out_vcf} ]; then
        echo "[$(date)] Skipping ${interval_name} — output already exists"
        return 0
    fi

    # check that the GenomicsDB for this interval exists
    if [ ! -d ${db_path} ]; then
        echo "[$(date)] ERROR: GenomicsDB not found for ${interval_name} at ${db_path}" >&2
        return 1
    fi

    echo "[$(date)] Running GenotypeGVCFs on interval ${interval_name}..."
    gatk GenotypeGVCFs \
        -R ${ref} \
        -V gendb://${db_path} \
        -O ${out_vcf} \
        -L ${interval_file} \
        --tmp-dir ${work_dir}/tmp \
        --include-non-variant-sites
        > ${log_file} 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date)] Finished GenotypeGVCFs for interval ${interval_name}"
    else
        echo "[$(date)] ERROR on GenotypeGVCFs for interval ${interval_name}" >&2
        return 1
    fi
}

export -f run_genotype
export work_dir ref db_dir genotype_dir log_dir

# ── 3. run GenotypeGVCFs in parallel across intervals
echo "[$(date)] Starting GenotypeGVCFs across all intervals..."

ls ${interval_dir}/*.interval_list | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog ${log_dir}/genotype_joblog.txt \
             --halt soon,fail=1 \
             run_genotype {}


# ── 4. GatherVcfs — merge all per-interval VCFs into a single final VCF
merged_vcf=${work_dir}/merged_vcf

mkdir -p ${merged_vcf}

final_vcf=${merged_vcf}/batch1_34samples.vcf.gz

# build sorted input list — must be in genome order
input_vcfs=$(ls ${genotype_dir}/*.genotyped.vcf.gz | sort -V | \
             awk '{print "--INPUT "$1}' | tr '\n' ' ')

gatk GatherVcfs \
    ${input_vcfs} \
    --OUTPUT ${final_vcf} \
    > ${log_dir}/gather_final.log 2>&1

# index the final VCF
gatk IndexFeatureFile \
    -I ${final_vcf} \
    >> ${log_dir}/gather_final.log 2>&1

if [ $? -eq 0 ]; then
    echo "[$(date)] Merging vcfs complete. Final VCF: ${final_vcf}"
else
    echo "[$(date)] ERROR in Merging vcfs — check ${log_dir}/gather_final.log" >&2
    exit 1
fi
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