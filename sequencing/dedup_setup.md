### Understanding how DeDup works and what it needs

DeDup needs to read Forward versus Reverse information in the header of the reads in the bams. 

In the second column of the bam file, there is a number. The value of this number actually tells us about the read orientation, 
if it is paired or not, pretty cool: https://broadinstitute.github.io/picard/explain-flags.html

In each of the bam, I thus need to create a new bam file that adds this information. 

To each of the merged bam read, I need to add M_

To each of the non-merged bam reads, I need to add R_ to reverse and F_ to forward strands. 

### Adding M_ to collapsed reads
```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

module load samtools
collapsed_bams=/scratch/midway3/rozennpineau/herbarium_partial_lane/collapsed_bams
out=/scratch/midway3/rozennpineau/herbarium_partial_lane/collapsed_bams/renamed_bams
threads=2

cd $collapsed_bams
for bam in *.sorted.bam; do

        #prefix header with M_
        name=${bam%.sorted.bam}
        samtools view -h ${name}.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/M_/' | samtools view -bS - > $out/${name}.prefixed.bam

        #sort and index
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o $out/${name}.prefixed.sorted.bam $out/${name}.prefixed.bam                         

done
```


For Forward and Reverse: (1) split between forward and reverse, then (2) add the prefix to the header. 


