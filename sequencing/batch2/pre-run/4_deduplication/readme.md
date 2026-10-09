## Prepping the reads for deduplication

### Prepping unmerged reads

script sent as a job array : first step is to make the list of files that need to go through the procedure

(1) generate list of files to prep

```
cd /scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged
for bam in *.unmerged.sorted.bam; do
echo $bam
done > unmerged_bams_to_rename.txt

```
(2) script name: prep_unmerged_F_R_bams_dedup.sh 

```
#!/bin/bash
#SBATCH --job-name=prep_bams_F_R
#SBATCH --output=logs/prep_bams_F_R_%a.out   # one log per sample
#SBATCH --error=logs/prep_bams_F_R_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G

# activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba
module load samtools

workdir=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged
threads=$SLURM_CPUS_PER_TASK

cd $workdir
mkdir -p forward reverse tmp logs

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# Each job picks one BAM from the todo list based on its array task ID.
# ─────────────────────────────────────────────────────────────────────────────
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $workdir/unmerged_bams_to_rename.txt)
name=${bam%.unmerged.sorted.bam}

echo "Array task $SLURM_ARRAY_TASK_ID processing: $name"

# ─────────────────────────────────────────────────────────────────────────────
# FORWARD READS
# -F 0x10 excludes reads with the reverse complement flag set,
# keeping only forward reads
# ─────────────────────────────────────────────────────────────────────────────
if samtools quickcheck forward/${name}.F.prefixed.unmerged.bam 2>/dev/null; then
    echo "forward/${name}.F.prefixed.unmerged.bam exists and is valid, skipping."
else
    echo "Splitting forward reads for $name..."
    samtools view -F 0x10 -h ${name}.unmerged.sorted.bam \
        | samtools view -bS - \
        > forward/${name}.F.unmerged.sorted.bam

    echo "Adding F_ prefix to read names for $name..."
    samtools view -h forward/${name}.F.unmerged.sorted.bam \
        | sed 's/^@/&/;/^[^@]/s/^/F_/' \
        | samtools view -bS - \
        > forward/${name}.F.prefixed.unmerged.bam
fi

echo "Sorting forward/${name}.F.prefixed.unmerged.bam..."
sambamba sort \
    -m 15GB \
    --tmpdir tmp \
    -t $threads \
    -o forward/${name}.F.prefixed.unmerged.sorted.bam \
    forward/${name}.F.prefixed.unmerged.bam

# ─────────────────────────────────────────────────────────────────────────────
# REVERSE READS
# -f 0x10 keeps only reads with the reverse complement flag set,
# keeping only reverse reads
# ─────────────────────────────────────────────────────────────────────────────
if samtools quickcheck reverse/${name}.R.prefixed.unmerged.bam 2>/dev/null; then
    echo "reverse/${name}.R.prefixed.unmerged.bam exists and is valid, skipping."
else
    echo "Splitting reverse reads for $name..."
    samtools view -f 0x10 -h ${name}.unmerged.sorted.bam \
        | samtools view -bS - \
        > reverse/${name}.R.unmerged.sorted.bam

    echo "Adding R_ prefix to read names for $name..."
    samtools view -h reverse/${name}.R.unmerged.sorted.bam \
        | sed 's/^@/&/;/^[^@]/s/^/R_/' \
        | samtools view -bS - \
        > reverse/${name}.R.prefixed.unmerged.bam
fi

echo "Sorting reverse/${name}.R.prefixed.unmerged.bam..."
sambamba sort \
    -m 15GB \
    --tmpdir tmp \
    -t $threads \
    -o reverse/${name}.R.prefixed.unmerged.sorted.bam \
    reverse/${name}.R.prefixed.unmerged.bam

echo "Finished $name."
```
(3) Submitting the job:
```
N=$(wc -l < /scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged/unmerged_bams_to_rename.txt)
echo "Submitting $N jobs"
sbatch --array=1-$N%8 prep_unmerged_F_R_bams.sh
```
(4) remove intermediate bams
This create a lot of intermediate files. 

The final files to keep are ".prefixed.unmerged.sorted.bam" and ".prefixed.unmerged.sorted.bam.bai". 
Remove "R.unmerged.sorted.bam" and "R.prefixed.unmerged.bam" when done. 

### Prepping collapsed reads

(1) generate list of files to prep

```
cd /scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed
for bam in *.collapsed.sorted.bam; do
echo $bam
done > collapsed_bams_to_rename.txt

```
(2) script name: prep_collapsed_M_bams.sh

```
#!/bin/bash
#SBATCH --job-name=prep_bams_M
#SBATCH --output=logs/prep_bams_M_%a.out   # one log per sample
#SBATCH --error=logs/prep_bams_M_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G

# activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba
module load samtools

workdir=/scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed
threads=$SLURM_CPUS_PER_TASK

cd $workdir
mkdir -p renamed_M tmp

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# Each job picks one BAM from the todo list based on its array task ID.
# ─────────────────────────────────────────────────────────────────────────────
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $workdir/collapsed_bams_to_rename.txt)
name=${bam%.collapsed.sorted.bam}

echo "Array task $SLURM_ARRAY_TASK_ID processing: $name"

if samtools quickcheck renamed_M/${name}.M.prefixed.bam 2>/dev/null; then
    echo "renamed_M/${name}.M.prefixed.bam exists and is valid, skipping."
else
    samtools view -h ${name}.collapsed.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/M_/' | samtools view -bS - > renamed_M/${name}.M.prefixed.bam
fi  

echo "Sorting renamed_M/${name}.M.prefixed.bam..."
sambamba sort \
    -m 15GB \
    --tmpdir tmp \
    -t $threads \
    -o renamed_M/${name}.M.prefixed.sorted.bam \
    renamed_M/${name}.M.prefixed.bam

echo "Finished $name."
```


(3) Submitting the job:
```
M=$(wc -l < /scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed/collapsed_bams_to_rename.txt)
echo "Submitting $M jobs"
sbatch --array=1-$M%8 prep_collapsed_M_bams.sh
```

After running the job once, 524 files were processes and I am missing 12. Identify those 12 samples: 

```
# Check each sample and keep only those that need processing
cd /scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed

for r1 in *.collapsed.sorted.bam; do
    prefx=${r1%.collapsed.sorted.bam}
    output=/scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed/renamed_M/${prefx}.M.prefixed.sorted.bam
    
    # If quickcheck fails (file missing or corrupt), add to todo list
    if ! samtools quickcheck $output 2>/dev/null; then
        echo $r1
    fi
done > collapsed_bams_to_rename.txt

```

Submit the job: 

```
N=$(wc -l < /scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed/collapsed_bams_to_rename.txt)
echo "Submitting $N jobs"
sbatch --array=1-$N%8 prep_collapsed_M_bams.sh
```

## Merge prefixed bams
(1) generate list of files to prep

```
cd /scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged/reverse
for bam in *.R.prefixed.unmerged.sorted.bam; do
echo $bam
done > bams_to_merge.txt #536 lines
```

```
#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --output=logs/merge_%a.out
#SBATCH --error=logs/merge_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G

# activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba
module load samtools

# set threads after environment is loaded
threads=$SLURM_CPUS_PER_TASK

unmerged_R_bams=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged/reverse
unmerged_F_bams=/scratch/midway2/rozennpineau/herbarium/batch2/bams/unmerged/forward
collapsed_bams=/scratch/midway2/rozennpineau/herbarium/batch2/bams/collapsed/renamed_M
output_folder=/scratch/midway2/rozennpineau/herbarium/batch2/bams/merged
tmp_dir=$output_folder/tmp

mkdir -p $output_folder $tmp_dir logs

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# ─────────────────────────────────────────────────────────────────────────────
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $unmerged_R_bams/bams_to_merge.txt)
name=${bam%.R.prefixed.unmerged.sorted.bam}

echo "Array task $SLURM_ARRAY_TASK_ID processing: $name"

# ─────────────────────────────────────────────────────────────────────────────
# INTEGRITY CHECK — skip if output already exists and is valid
# ─────────────────────────────────────────────────────────────────────────────
if samtools quickcheck $output_folder/${name}.prefixed.sorted.bam 2>/dev/null; then
    echo "$name already merged and sorted, skipping."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# INPUT CHECK — verify all three input BAMs exist before merging
# Missing inputs would cause samtools merge to fail silently or partially
# ─────────────────────────────────────────────────────────────────────────────
R_bam=$unmerged_R_bams/${name}.R.prefixed.unmerged.sorted.bam
F_bam=$unmerged_F_bams/${name}.F.prefixed.unmerged.sorted.bam
M_bam=$collapsed_bams/${name}.M.prefixed.sorted.bam

for input in "$R_bam" "$F_bam" "$M_bam"; do
    if ! samtools quickcheck $input 2>/dev/null; then
        echo "ERROR: missing or corrupted input: $input"
        exit 1
    fi
done

# ─────────────────────────────────────────────────────────────────────────────
# MERGE
# ─────────────────────────────────────────────────────────────────────────────
echo "Merging BAMs for $name..."
samtools merge \
    $output_folder/${name}.prefixed.bam \
    $R_bam \
    $F_bam \
    $M_bam

# ─────────────────────────────────────────────────────────────────────────────
# SORT
# ─────────────────────────────────────────────────────────────────────────────
echo "Sorting merged BAM for $name..."
sambamba sort \
    -m 15GB \
    --tmpdir $tmp_dir \
    -t $threads \
    -o $output_folder/${name}.prefixed.sorted.bam \
    $output_folder/${name}.prefixed.bam

echo "Finished $name."

```

Submit the job:
```
n=$(wc -l < $unmerged_R_bams/bams_to_merge.txt)
echo "Submitting $n jobs"
sbatch --array=1-$n%8 merge_bams_job_array.sh
```

## Split into scaffolds

(1) Prepare list of bams to split
```

cd /scratch/midway3/rozennpineau/herbarium/batch2/bams/merged
for bam in *.prefixed.sorted.bam; do
    echo $bam
done > bams_to_split.txt #536 lines

```
(2) Prepare script
```
#!/bin/bash
#SBATCH --job-name=split_scaffolds
#SBATCH --output=logs/split_%a.out
#SBATCH --error=logs/split_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G

# activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
module load samtools

bams=/scratch/midway3/rozennpineau/herbarium/batch2/bams/merged
workdir=/scratch/midway2/rozennpineau/herbarium/batch2/bams/scaffolded

mkdir -p $workdir logs
cd $workdir

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# ─────────────────────────────────────────────────────────────────────────────
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $bams/bams_to_split.txt)
name=${bam%.prefixed.sorted.bam}

echo "Array task $SLURM_ARRAY_TASK_ID processing: $name"

# ─────────────────────────────────────────────────────────────────────────────
# INPUT CHECK
# ─────────────────────────────────────────────────────────────────────────────
if ! samtools quickcheck $bams/$bam 2>/dev/null; then
    echo "ERROR: missing or corrupted input: $bams/$bam"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# COMPLETION CHECK
# ─────────────────────────────────────────────────────────────────────────────
flag=$workdir/$name/split.done

if [ -f "$flag" ]; then
    echo "$name has already been split, skipping."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# SPLIT INTO SCAFFOLDS
# samtools idxstats returns one line per scaffold plus a final '*' wildcard
# line for unmapped reads — head -n -1 drops that last line.
# Each scaffold is written to its own BAM and indexed.
# ─────────────────────────────────────────────────────────────────────────────
mkdir -p $workdir/$name

for chr in $(samtools idxstats $bams/$bam | cut -f1 | head -n -1); do
    echo "Splitting $chr from $name..."
    samtools view -b $bams/$bam "$chr" > $workdir/$name/${chr}.bam
    samtools index $workdir/$name/${chr}.bam
done

# ─────────────────────────────────────────────────────────────────────────────
# WRITE FLAG FILE
# Only reached if the loop above completed without error.
# Acts as a reliable completion marker for future runs.
# ─────────────────────────────────────────────────────────────────────────────
echo "Split completed on $(date)" > $flag
echo "Finished splitting $name."

```
(3) Submit job
```
n=$(wc -l < $bams/bams_to_split.txt)
echo "Submitting $n jobs"
sbatch --array=1-$n%8 split_bams_to_scaffolds_job_array.sh
```

This should create a folder per sample, with 1060 bams in each folder for each scaffold.

## Run Dedup

(1) Prepare list of bams to split
```
cd /scratch/midway2/rozennpineau/herbarium/batch2/bams/scaffolded
for sample in sample_*; do
    echo $sample
done > samples_to_dedup.txt #536 lines
```
(2) Prepare script

```
#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --output=logs/dedup_%a.out
#SBATCH --error=logs/dedup_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2GB

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/dedup/

module load samtools

work_dir=/scratch/midway2/rozennpineau/herbarium/batch2/bams/scaffolded
cd $work_dir

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# ─────────────────────────────────────────────────────────────────────────────
samp=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $work_dir/samples_to_dedup.txt)
cd $samp
for bam in *.bam; do
    name=${bam%.bam}
    mkdir -p $name
    dedup -i $bam -o ${name}
done
```
(3) Submit job
```
sbatch --array=1-536%8 run_dedup_job_array.sh
```

## Merge bams after DeDup

```
start_dir=/scratch/midway2/rozennpineau/herbarium/batch2/bams/scaffolded

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bamtools
module load samtools

ulimit -n 4096 # increase upper limit of number of files that can be opened at once

cd $start_dir
for dir in ./*; do #list directories one level down only 
    cd $dir

    realpath */[Ss]*rmdup.bam > bams_to_merge.list
    echo "merging $dir..."
    bamtools merge -list bams_to_merge.list -out $dir.scaffolds.dedup.bam
    samtools index $dir.scaffolds.dedup.bam #index

    cd ..
done
```

## Calculate Duplication rate

```
echo -e "Sample\tScaffold1-16_duplication_rate" > scaffold1-16_mean_dup_rate.txt

for dir in ./sample_*; do
  #echo $dir
  dup_rate=$(cat $dir/dup_rate_summary.txt | grep Scaffold | awk 'NR>1 {sum += $2; n++} END {print sum/n}')
  echo -e "${dir}\t${dup_rate}" >> Scaffold1-16_mean_dup_rate.txt
done
```

duplication rates for samples on the one hand, and number of bp per mapped read before dedup on the other
send those to Julia

## Calculate coverage after deduplication

Run as a job array. 

```
#!/bin/bash
#SBATCH --job-name=calc_bp
#SBATCH --output=logs/calc_bp_%a.out
#SBATCH --error=logs/calc_bp_%a.err
#SBATCH --time=02:00:00              
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

module load samtools

path=/scratch/midway2/rozennpineau/herbarium/batch2/bams/dedup
mkdir -p $path/logs $path/coverage_results
cd $path

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE SELECTION
# ─────────────────────────────────────────────────────────────────────────────
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $path/bams_to_calc.txt)
name=${bam%.scaffolds.dedup.bam}

echo "Array task $SLURM_ARRAY_TASK_ID processing: $name"

if ! samtools quickcheck $bam 2>/dev/null; then
    echo "ERROR: missing or corrupted input: $bam"
    exit 1
fi

# # ─────────────────────────────────────────────────────────────────────────────
# CALCULATE TOTAL BASE PAIRS MAPPED
# "bases mapped (cigar)" is the most accurate measure of base pairs mapped —
# it counts only bases that are part of the actual alignment, excluding
# soft-clipped bases that are in the read but not aligned to the reference.
# ─────────────────────────────────────────────────────────────────────────────
bases=$(samtools stats -@ $SLURM_CPUS_PER_TASK $bam \
    | grep "^SN" \
    | grep "bases mapped (cigar):" \
    | cut -f 3)

echo -e "$name\t$bases" > $path/coverage_results/${name}.bp.txt
echo "Finished $name: $bases base pairs mapped."

```

Starting job :

```
sbatch --array=1-536%8 calc_bp_after_dedup.sh

```

Combine all samples in one file: 

```
cat $path/coverage_results/*.coverage.txt >> nb_bp_after_dedup.txt
```