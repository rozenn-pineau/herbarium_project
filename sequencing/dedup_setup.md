### Understanding how DeDup works and what it needs

DeDup needs to read Forward versus Reverse information in the header of the reads in the bams. 

In the second column of the bam file, there is a number. The value of this number actually tells us about the read orientation, 
if it is paired or not, pretty cool: https://broadinstitute.github.io/picard/explain-flags.html

In each of the bam, I thus need to create a new bam file that adds this information. 

To each of the merged bam read, I need to add M_

To each of the non-merged bam reads, I need to add R_ to reverse and F_ to forward strands. 

