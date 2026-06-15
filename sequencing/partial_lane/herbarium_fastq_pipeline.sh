# Software requirements:
#	AdapterRemoval: https://github.com/mikkelschubert/adapterremoval
#	bwa: https://github.com/lh3/bwa
#	samtools: https://www.htslib.org/download/
#	picard: https://github.com/broadinstitute/picard/releases/
# MapDamage2: https://ginolhac.github.io/mapDamage/
 #
# USAGE:
# bash script.sh <reference_genome.fasta> <n_threads>
#
# You can redirect all stderr to /dev/null if needed. bash script.sh ...options... 2> /dev/null

# Change paths accordingly
#bwa=/path/to/bwa
#samtools=/path/to/samtools
fastp=/ohta/julia.kreiner/software/fastp
picard=/ohta/julia.kreiner/software/picard.jar
#path=/ohta/julia.kreiner/waterhemp/herbarium
prefx=$3
pathtobams=/ohta2/julia.kreiner/waterhemp/herbarium/fastqs

reference=$1
threads=$2

# Assuming names are indivexmple1_R1.fastq.gz, indivexmple1.R2.fastq.gz
# It keeps whatever is before the first dot '.' as basename for the downstream outputs. You can change if needed

cd $path/raw_data/$prefx

cat ${prefx}*_1.fq.gz > ${prefx}.R1.fastq.gz
cat ${prefx}*_2.fq.gz > ${prefx}.R2.fastq.gz

# 1. Remove adapters, polyQ tails, and merge reads (important for short frags)
$fastp --in1 ${prefx}.R1.fastq.gz --in2 ${prefx}.R2.fastq.gz --out1 ${prefx}.R1.unmerged.fastq.gz --out2 ${prefx}.R2.unmerged.fastq.gz  --merge --merged_out ${prefx}.collapsed.gz

# 2. Map merged (collapsed) reads to Reference Genome and calculate Endougenous DNA
bwa mem -t $threads -R "@RG\tID:$prefx\tSM:$prefx" $reference ${pathtobams}/${prefx}_R1.fastq.gz ${pathtobams}/${prefx}_R2.fastq.gz | samtools view -@ $threads -Sbh - >  /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/${prefx}.uns.bam

#cd $path/bams

total=$(samtools view -@ $threads -c $prefx.full.uns.bam)
echo -e "TotalReads\n$total" >> $prefx.log

sambamba sort -m 15GB --tmpdir $path/bams/tmp -t $threads -o /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/$prefx.sorted.scaled.bam /ohta2/julia.kreiner/waterhemp/herbarium/femaleref/${prefx}.scaled.bam
samtools merge ${prefx}.final.sorted.bam ${prefx}.unmerg.sorted.bam ${prefx}.merged.sorted.bam
rm ${prefx}.uns.bam
mapped=$(samtools view -@ $threads -c $prefx.bam)
echo -e "MappedReads\n$mapped" >> ${prefx}.log
echo "EndogenousDNA" >> ${prefx}.log
python -c "print(float($mapped)/ $total)" >> ${prefx}.log

# 2b. MapDamage to calculate deamination levels per sample, across the reads, and rescale per-base quality
mapDamage -i ${prefx}.final.sorted.bam -r myreference.fasta --rescale


#3. Mark duplicates with DeDup - optimized for merged pair-end reads (aDNA)

mkdir {}; java -jar ~/software/DeDup/DeDup-0.12.6.jar -i /ohta2/julia.kreiner/herbarium/femaleref/bams/{prefx}_rescaled/{}.final.sorted.bam -m -o {}"



# 4. Remove non-marked file and generate index
mv $prefix.dd.bam $prefix.bam
samtools index /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.dd.bam
samtools flagstat /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.dd.bam > /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.stats

java -jar $picard ValidateSamFile I=/ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.final.sorted_rmdup.bam MODE=SUMMARY O=/ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.check
#
