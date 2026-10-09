## Using Gatk to call variants 

Overview of steps for variant calling with GATK:

1. Split the reference genome into windows. 
In our case, the genome is 650 MB and here I have 108 files. If I split in 100 KB windows, that 6500 shards. If I split in 5 MB windows, that 130 shards. I am trying to split in 5 MB windows for now to be able to split the job as job arrays and have the node work on several jobs as a time with GNU parallel as well. For 5 MB windows, this is already 130 * 108 sub nodes working.  

GATK is written in Java, so every gatk command starts a Java Virtual Machine, which then runs GATK’s code inside it. Each GATK launch has a startup cost (loading the JVM and the tool), hence running many short GATK jobs is wasteful.


2. Run HaplotypeCaller independently on each BAM to produce a gVCF.
A gVCF has information about all loci in the samples, not just the variant ones. You can choose to use an option to have every single locus in a single line, with the read information for each locus (--ERC BP_RESOLUTION) or to have blocks of similar coverage summarized in blocks (--ERC GVCF).

3. Combine the gVCFs.

4. Run GenomicsDBImport. 
GenomicsGDImport turns the sample-centered file into a genome-centered file: before there is one file per sample. Now there is one file per genomic region that combines all samples. 

- Joint genotype all samples.

### Step 0. Getting the bams for the 2022 herbarium dataset
Logon to the NAS, then transfer bams to Midway 3. 

```
scp -r *bam* rozennpineau@midway3.rcc.uchicago.edu:/scratch/midway3/rozennpineau/herbarium/2022/bams/

```
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


Reference and dictionary are here: /project/kreiner/data/genome/ 

### Step 2. Split the reference genome into 5MB intervals 

A shard is one slice of the genome that gets processed as a separate job. Initially, I tried splitting the genome in 100 kb windows. However, the script launches one HaplotypeCaller per BAM per shard. With ~6,547 shards of 100 kb, that’s:

110 × 6,547 ≈ 720,000 gatk launches

Each launch pays JVM startup and reference loading, so most of the runtime would be overhead. to fix this, I am trying a bigger shar of 5 Mb (--scatter-count 130) for HaplotypeCaller. That’s a reasonable balance: each run lasts long enough to amortize the startup cost, and you can still resume at a fine granularity.

Updated code: 
```
#!/bin/bash
#SBATCH --job-name=make_gvcfs
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --time=36:00:00
#SBATCH --array=1-1%40     # one task per BAM, 40 running at once
#SBATCH --output=logs/make_gvcfs_%A_%a.out
#SBATCH --error=logs/make_gvcfs_%A_%a.err

module load gatk samtools parallel

#load directories
main_dir=/scratch/midway3/rozennpineau/herbarium/2022
bam_dir=${main_dir}/bams
work_dir=${main_dir}/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
log_dir=${work_dir}/logs
mkdir -p ${work_dir}/tmp ${log_dir}

# one sample per array task
bam=$(ls ${bam_dir}/*_rmdup_rescaled_resorted.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(basename ${bam} _rmdup_rescaled_resorted.bam)
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

echo "[$(date)] ${sample} complete."
```



### Step 3. HaplotypeCaller - make GVCF files

!!! Next time, add the merging gvcfs step at the end of this file. 
-- no actually, keep steps separate to be able to check them. An error occuring in the gvcf making might be invisible at the gvcf merging step.

This is using parallel to run more intervals at a time. I could transform this code to run in as a job array, one interval at a time. 
Check how much time it took for one node, 9 jobs per node (each sample is divided in 100kb windows).


script name: run_gatk_1.sh

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
#SBATCH --array=2-110%40     # one task per BAM, 40 running at once
#SBATCH --output=logs/make_gvcfs_%A_%a.out
#SBATCH --error=logs/make_gvcfs_%A_%a.err

module load gatk samtools parallel

#load directories
main_dir=/scratch/midway3/rozennpineau/herbarium/2022
bam_dir=${main_dir}/bams
work_dir=${main_dir}/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
log_dir=${work_dir}/logs
mkdir -p ${work_dir}/tmp ${log_dir}

# one sample per array task
bam=$(ls ${bam_dir}/*_rmdup_rescaled_resorted.bam | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(basename ${bam} _rmdup_rescaled_resorted.bam)
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
I tested the script on one bam and checked the memory allocation and job efficiency:
```
sacct -j 60323231 --format=JobID,State,Elapsed,MaxRSS,AllocCPUS,TotalCPU
```
The job took 2h16 to complete for one bam. I also used 14GB.  
With the current code, I start 40 jobs at a time and I decreased the memory allocation to 24GB. 

on the remaining 107 BAMs: job 60333290
Check that all 107 jobs were completed:
```
sacct -j 60333290 -X -n --format=State | grep -c COMPLETED     # expect 107
```
### Check for missing shard/gvcfs
```
work_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
n_expected=$(ls ${interval_dir}/*.interval_list | wc -l)

i=0
for bam in ${bam_dir}/*_rmdup_rescaled_resorted.bam; do
    i=$((i+1))
    sample=$(basename ${bam} _rmdup_rescaled_resorted.bam)
    n_ok=$(ls ${work_dir}/${sample}/*.g.vcf.gz.tbi 2>/dev/null | wc -l)
    [ "$n_ok" -eq "$n_expected" ] || echo "index ${i}  ${sample}: ${n_ok}/${n_expected}"
done
```
Looking good !

### Step 4 - Merging gvcf files

```
# modules
module load gatk
module load parallel

# ── 1. define variables
work_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf
unmerged_dir=${work_dir}/unmerged_gvcfs
merged_dir=${work_dir}/merged_gvcfs
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
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
ls -d ${unmerged_dir}/HB*/ | \
    parallel --jobs ${SLURM_CPUS_PER_TASK} \
             --joblog ${log_dir}/gather_joblog.txt \
             --halt soon,fail=1 \
             gather_sample {}
```

Check if resources were efficiently used: 60346610

```
sacct -j 60312659 --format=JobID,JobName,User,AllocCPUS,ReqMem,ReqTRES,Elapsed,MaxRSS
```

Thi shows that 2000K were used and I requested 10GB of memory so I need much much less. I will change the script to reflect that.

### Step 5 - GenomicsDBImport

What does [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) do: 
Joint genotyping with GenotypeGVCFs needs to consider all samples simultaneously, but reading a lot (500) GVCFs directly is slow and inefficient. GenomicsDBImport solves this by consolidating them into an optimized database format. It reorganizes the data from being sample-oriented to site-oriented, to then be able to call variantas across all samples at each position in the genome.


Options for GenomicsDBImport that we night want/need to play with:
--batch-size 0	Batch size controls the number of samples for which readers are open at once and therefore provides a way to minimize memory consumption. However, it can take longer to complete. Use the consolidate flag if more than a hundred batches were used. This will improve feature read time. batchSize=0 means no batching (i.e. readers for all samples will be opened at once) Defaults to 0.

--reader-threads
1	How many simultaneous threads to use when opening VCFs in batches; higher values may improve performance when network latency is an issue


1. Create the sample map:

```
work_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf
merged_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf/merged_gvcfs
sample_map=${work_dir}/sample_map.txt
> ${sample_map}
for gvcf in ${merged_dir}/*.merged.g.vcf.gz; do
    echo -e "$(basename ${gvcf} .merged.g.vcf.gz)\t${gvcf}" >> ${sample_map}
done
wc -l ${sample_map}      # expect 108
```

2. Run GenomicsDBImport

```
#!/bin/bash
#SBATCH --job-name=gdbimport
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00           
#SBATCH --array=1-130%50
#SBATCH --output=logs/gdb_%A_%a.out
#SBATCH --error=logs/gdb_%A_%a.err

module load gatk

# load directories
work_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf
merged_dir=${work_dir}/merged_gvcfs
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB
db_dir=${work_dir}/genomicsdb
log_dir=${work_dir}/logs
sample_map=${work_dir}/sample_map.txt

mkdir -p ${db_dir} ${log_dir} ${work_dir}/tmp

# check that the sample map has been created
if [ ! -s ${sample_map} ]; then
    echo "ERROR: ${sample_map} missing; create it before submitting" >&2
    exit 1
fi

# one shard per array task
interval_file=$(ls ${interval_dir}/*.interval_list | sed -n "${SLURM_ARRAY_TASK_ID}p") # pick the Nth shard file for array task N.
interval_name=$(basename ${interval_file} .interval_list)
db_path=${db_dir}/${interval_name}
tmp_dir=${work_dir}/tmp/${interval_name}
mkdir -p ${tmp_dir}

# skip if finished; clean up if a previous attempt was interrupted
if [ -f ${db_path}/callset.json ]; then
    echo "${interval_name} already done"; exit 0
fi
rm -rf ${db_path}

gatk --java-options "-Xmx8g" GenomicsDBImport \
    --sample-name-map ${sample_map} \
    --genomicsdb-workspace-path ${db_path} \
    -L ${interval_file} \
    --merge-input-intervals true \
    --batch-size 50 \
    --reader-threads 2 \
    --tmp-dir ${tmp_dir} \
    > ${log_dir}/${interval_name}.genomicsdb.log 2>&1

status=$?
rm -rf ${tmp_dir}
exit ${status}

# Explanations of the command line choices:
# -Xmx sets the maximum heap size. -Xmx8g means “this JVM may use up to 8 GB of heap” - heap = memory

--merge-input-intervals true	Treats the intervals in the shard as one combined region and stores it as a single array, instead of creating one array per interval. This matters because the shards contain many small contigs.
--batch-size 50	The number of samples whose GVCFs are opened and read at the same time.
--reader-threads 2	Number of threads used to read GVCF files. More threads can speed things up, especially with slow storage. I should match --cpus-per-task.

```
Running on a single job first to adjust memory allocations: 60347356


```
sacct -j 60347356 --format=JobID,State,Elapsed,MaxRSS,AllocCPUS,TotalCPU
```
This showed that I could lower to memory to 8Gb and one task per node. I modified the header accordingly. 


Running on 130 shards: 60349280

The full job took about 2-3h (running 50 jobs at a time).
Check for any task that didn't complete cleanly:
```
sacct -j 60349280 -X --format=JobID,State,Elapsed,MaxRSS | grep -v COMPLETED #everything worked
```
Looking good so far here ! 

I will start genotyoe calling for this set of file to test the script. 

### Step 6 - Genotype GVCFs

[GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/21905118377755-GenotypeGVCFs)

```
#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --account=pi-kreiner
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=14G                 
#SBATCH --time=6:00:00           
#SBATCH --array=2-130%50
#SBATCH --output=logs/geno_%A_%a.out
#SBATCH --error=logs/geno_%A_%a.err

module load gatk

#load directories
work_dir=/scratch/midway3/rozennpineau/herbarium/2022/gvcf
ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
interval_dir=/scratch/midway3/rozennpineau/herbarium/ref_intervals_5MB   # same shards as the import step
db_dir=${work_dir}/genomicsdb
genotype_dir=${work_dir}/genotypes
log_dir=${work_dir}/logs

mkdir -p ${genotype_dir} ${log_dir} ${work_dir}/tmp

# one shard per array task
interval_file=$(ls ${interval_dir}/*.interval_list | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "${interval_file}" ]; then
    echo "ERROR: no interval file for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1
fi
interval_name=$(basename ${interval_file} .interval_list)
db_path=${db_dir}/${interval_name}
out_vcf=${genotype_dir}/${interval_name}.genotyped.vcf.gz
tmp_dir=${work_dir}/tmp/geno_${interval_name}

# skip if finished (index is written last)
if [ -f ${out_vcf}.tbi ]; then
    echo "${interval_name} already done"; exit 0
fi

# the database must be complete
if [ ! -f ${db_path}/callset.json ]; then
    echo "ERROR: no finished GenomicsDB for ${interval_name}" >&2; exit 1
fi

rm -f ${out_vcf} ${out_vcf}.tbi     # leftovers from an interrupted run
mkdir -p ${tmp_dir}

gatk --java-options "-Xmx8g" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${db_path} \
    -O ${out_vcf} \
    -L ${interval_file} \
    --tmp-dir ${tmp_dir} \
    --include-non-variant-sites \
    --max-alternate-alleles 4 \
    > ${log_dir}/${interval_name}.genotype.log 2>&1

status=$?
rm -rf ${tmp_dir}
exit ${status}

#Explanations
#-V gendb:// prefix tells GATK to read from a GenomicsDB workspace instead of a VCF file

```
```
sacct -j 60354943 -X --format=JobID,JobName,User,AllocCPUS,ReqMem,ReqTRES,Elapsed,MaxRSS
```

