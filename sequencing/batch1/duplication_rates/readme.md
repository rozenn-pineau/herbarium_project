### Removing duplicates

We do so using DeDup. 

DeDup needs to read Forward versus Reverse information in the header of the reads in the bams.
To each of the merged bam read, I need to add M_

To each of the non-merged bam reads, I need to add R_ to reverse and F_ to forward strands.

### Adding M_ to merged (collapsed reads)
```
### Adding M_ to merged (collapsed_) reads

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba
module load samtools

collapsed_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/collapsed
out=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/collapsed/renamed

threads=2

cd $collapsed_bams

for bam in *.sorted.bam; do

        #prefix header with M_
        name=${bam%.sorted.bam}
        samtools view -h ${name}.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/M_/' | samtools view -bS - > $out/${name}.prefixed.bam

        #sort and index
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o $out/${name}.prefixed.sorted.bam $out/${name}.prefixed.bam                         

done



### Adding F_ to forward reads
cd /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged
mkdir forward

for bam in *.unmerged.sorted.bam; do
	name=${bam%.unmerged.sorted.bam}
	# split
	samtools view -F 0x10 -h ${name}.unmerged.sorted.bam | samtools view -bS - > forward/${name}.F.unmerged.sorted.bam #exclude reverse read
	# rename
	samtools view -h forward/${name}.F.unmerged.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/F_/' | samtools view -bS - > forward/${name}.F.prefixed.unmerged.bam
	#sort and index
	sambamba sort -m 15GB --tmpdir tmp -t $threads -o forward/${name}.F.prefixed.unmerged.sorted.bam forward/${name}.F.prefixed.unmerged.bam   

done


### Adding R_ to reverse reads
cd /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged
mkdir reverse

for bam in *.unmerged.sorted.bam; do
        name=${bam%.unmerged.sorted.bam}
        # split
        samtools view -f 0x10 -h ${name}.unmerged.sorted.bam | samtools view -bS - > reverse/${name}.R.unmerged.sorted.bam #exclude forward read
        # rename
        samtools view -h reverse/${name}.R.unmerged.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/R_/' | samtools view -bS - > reverse/${name}.R.prefixed.unmerged.bam
        #sort and index
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o reverse/${name}.R.prefixed.unmerged.sorted.bam reverse/${name}.R.prefixed.unmerged.bam                                                          

done



```













