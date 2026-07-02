Align reads to reference genome using baw mem. 

Align unmerged and collapsed reads separately. 

### Aligning unmerged reads
```
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bwa
module load samtools

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams
threads=6
working_dir=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc_twice/sequenced

cd $working_dir


for r1 in *_R1.unmerged.fq.gz; do
        prefx=${r1%_R1.unmerged.fq.gz}
        output=$out/${prefx}.unmerged.uns.bam

        if [ -s "$output" ]; then #if output already exists and is greater than 0 bytes
                echo "$output already exists, skipping."
        else
                # Map unmerged reads to Reference Genome, turn SAM to BAM
                bwa mem -t $threads -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" $ref ${prefx}_R1.unmerged.fq.gz ${prefx}_R2.unmerged.fq.gz | samtools view -@ $threads -Sbh - >  $output
        fi
done


#explanation: align (produces sam files) bwa mem -t (number of threads) -R (read group header line) $reference_genome reads_lane1 reads_lane2
#pipe: converts sam to bams with samtoolis view @ #threads -Sbh sam file input, convert to bam, with header
```
I added a condition that checks for the presence of the file before running - this is because with 36 hours of wall time, about half of the samples had been aligned and I needed to start the script again. 

### Align collapsed reads

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bwa
module load samtools

ref=/project/kreiner/data/genome/Atub_193_hap2.fasta
out=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams
threads=6
working_dir=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc_twice/sequenced

cd $working_dir

for r1 in *.collapsed.fq.gz; do
        prefx=${r1%.collapsed.fq.gz}
        output=$out/${prefx}.uns.bam

        if [ -s "$output" ]; then #if output already exists and is greater than 0 bytes
                echo "$output already exists, skipping."
        else
        # Map unmerged reads to Reference Genome, turn SAM to BAM
        bwa mem -t $threads -R "@RG\tID:${prefx}\tSM:${prefx}\tPL:ILLUMINA\tLB:${prefx}" $ref ${prefx}.collapsed.fq.gz | samtools view -@ $threads -Sbh - >  $output
        fi
done

```





