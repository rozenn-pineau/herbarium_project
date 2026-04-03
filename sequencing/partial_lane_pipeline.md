#Goal: to record the pipeline for the partial lane sequencing. Draw from Julia's GitHub + Science paper.

Steps should include:
Comparison of C to T deamination between Julia's raw reads and the new raw reads using MapDamage. 
  --> 	does the lib prep decrease deamination?

Comparison before and after merging reads
  --> map before and after merging reads: set of reads paired and merged, map twice, merge resulting bams, compare results

Comparison of the number of duplicates between sets of herbarium sequence data AND between sample quality/concentrations. 

Compare mean fragment size between between sets of herbarium sequence data AND between sample quality/concentrations. 

tool: fastp
use: remove adapters, poly Q tails, merge reads
command line: 

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



# 4. Remove non-marked file and generate index
mv $prefix.dd.bam $prefix.bam
samtools index /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.dd.bam
samtools flagstat /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.dd.bam > /ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.stats

java -jar $picard ValidateSamFile I=/ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.final.sorted_rmdup.bam MODE=SUMMARY O=/ohta2/julia.kreiner/herbarium/femaleref/bams/${prefx}.check
#
