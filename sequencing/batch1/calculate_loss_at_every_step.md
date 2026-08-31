
I want tot know how much data is lost at every step to be able to calculate how many samples we can have per lane for the next rounds of sequencing. 


One challenge is that the coverage is sometimes in file size (Gbytes) and other times in number of bases (Gbases).
I am going back tp the raw files to caculate the number of bases in the fastq files (raw and after demultiplexing).

To do so, I use a tool named SeqKit.

```
conda activate /project/kreiner/rpineau/seqkit
conda install bioconda::seqkit

```

Script to use seqkit on the 125 fasq files:

```

```