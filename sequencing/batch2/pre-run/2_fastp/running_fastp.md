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
#SBATCH --nodes=4
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