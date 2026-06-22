### Doanloading the data from Novogene usinf lftp
```
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC26058003-Z01-F002:uD2dCtY9@usftp22.novogene.com; mirror --verbose --use-pget-n=8 -c'
```
### Checking the download 
```
md5sum -c MD5.txt
```

Every file looked good.

### Demultiplexing the fastq

(a) Match barcode to sequence:

I used a custom R script to match the barcodes i5 and i7 to their corresponding sequence (match_barcode_sequence.R).

(b) demultiplex:

I am choosing not to use [Ultraplex](https://github.com/ulelab/ultraplex) because it removes "bad quality bases" from the 3' end.
I tried flexbar () but it did not handle the two barcode systems very well. 
I am now trying [demultiplex](https://github.com/EichlerLab/demultiplex).

Installing demultiplex:
```
conda create -p /project/kreiner/rpineau/demultiplex
conda activate /project/kreiner/rpineau/demultiplex
git clone https://github.com/jfjlaros/demultiplex.git
cd demultiplex

pip install .

```

Using demultiplex command: 

usage: demultiplex demux [-h] [-r] [--format {normal,x,umi,unknown}] [-s START] [-e END] [-m MISMATCH] [-d]
                         [-p PATH]
                         BARCODES INPUT [INPUT ...]


I will try with the command line usually for long read data that aligns the read to look for the barcode:
```

conda activate /project/kreiner/rpineau/demultiplex

#working on smaller files
zcat P1_WKDL260007250-1A_23C5CJLT4_L4_1.fq.gz | head -n 80 > head_L1.fq.gz
zcat P1_WKDL260007250-1A_23C5CJLT4_L4_2.fq.gz | head -n 80 > head_L2.fq.gz 

#demultiplex match command
demultiplex match output_barcode_sequences.txt head_L1.fq head_L2.fq
# error message: AttributeError: 'str' object has no attribute '_accessor'
```


I could not fix the tool's error so I am trying a new approach with a script I found online. 

Testing on smaller files:
```
#activate java
conda activate /home/rozennpineau/java
conda activate /home/rozennpineau/bbmap

R1=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/head_L1.fq.gz
R2=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/head_L2.fq.gz
barcodes=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxbyname_batch1_barcodes.txt
out_folder=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed

./demuxbyname.sh in=$R1 in2=$R2 \
  out=$out_folder/sample_%_R1.fastq.gz out2=$out_folder/sample_%_R2.fastq.gz \
  delimiter=: prefixmode=f \
  names=$barcodes hdist=1 \
  outu=$out_folder/unmatched_R1.fastq.gz outu2=$out_folder/unmatched_R2.fastq.gz \
  -Xmx16g

# Trying on larger files:
zcat P1_WKDL260007250-1A_23C5CJLT4_L4_1.fq.gz | head -400000 | gzip > head_L1.fq.gz 
zcat P1_WKDL260007250-1A_23C5CJLT4_L4_2.fq.gz | head -400000 | gzip > head_L2.fq.gz 

```
From this test run, a lot of barcodes were sent to the unmatch pile - I will try a first run with Hamming distance of 1. I may need to change tool or hamming distance to 2 to increase yield.


This is not working really well. I am trying to find a solution that leverages the barcodes from the reads themselves. 
Some checks.
Where are the barcodes in the sequence?
```
#20 instances of i7 in R1 reads, 0 in R2 reads
grep -o -E 'ACGTTACC' first1000_R1.txt | sort | uniq -c
#1 instance of this i5 in R1 and R2 reads
grep -o -E 'ACCGACAA' first1000_R2.txt | sort | uniq -c 

grep -n 'AGTGGCAA' first1000_R2.txt
214:TAAACGTTTTTTTTAAGAAAGTGTGCATTTGTGTGCACAAAACAAAAAAGATAACGTTTTCAAAATTTCTACTTGATGCAGTGAATTGATGAATGCAATAATCTTTCCATGAGTGGCAATAGGGAAAAACTATTCGGTTTGAGGTGTTCT
495:GTAAGGCCGAGCATGTCAGTGGCAAGCTATCTATTCCTAGTATTAAAAAACAACATATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGTCTGTGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAGGGGGGGGGGGGGGGG
726:ATTGCCACATTTATCAACTCGAGCATGAACCTTTTTGGAGCTTCTTTCTTAGGTATGAGACGTTCTTCCTCCGCATACGTAATGTTGAAAAGTGGCAAATTTGTCACCATCTTTATGGAGGATTGAATGTTGAGACTAGGGAACTCGTTG

grep -n 'AGTGGCAA' first1000_R2.txt
214:TAAACGTTTTTTTTAAGAAAGTGTGCATTTGTGTGCACAAAACAAAAAAGATAACGTTTTCAAAATTTCTACTTGATGCAGTGAATTGATGAATGCAATAATCTTTCCATGAGTGGCAATAGGGAAAAACTATTCGGTTTGAGGTGTTCT
495:GTAAGGCCGAGCATGTCAGTGGCAAGCTATCTATTCCTAGTATTAAAAAACAACATATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGTCTGTGGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAGGGGGGGGGGGGGGGG
726:ATTGCCACATTTATCAACTCGAGCATGAACCTTTTTGGAGCTTCTTTCTTAGGTATGAGACGTTCTTCCTCCGCATACGTAATGTTGAAAAGTGGCAAATTTGTCACCATCTTTATGGAGGATTGAATGTTGAGACTAGGGAACTCGTTG
(base) [rozennpineau@midway3-login6 P1]$ grep -n 'ACGTTACC' first1000_R2.txt 
(base) [rozennpineau@midway3-login6 P1]$ grep -n 'ACGTTACC' first1000_R1.txt 
160:TNCTTTTTAACATTTTTTGTTTTAGACATAAAATTCAAGTCGTTTTGCTTTGTACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTAGGGTGGCATCTTCTGCTTTTAAAAAGGGGGGGGGGGGGGGGGGGGGGG
171:CNCTATGAAACTGTGGGAAAGAGTTATTGAGAGAAGAGTTAGACGGGAGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTAGTATGTCATCTTCTTCTTTCAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGG
195:CNAACTCAAATCTGACCGCTCTCTCCATACGAACTCCGAAAGATGATATGATGAGAAGGCCTGCCCCACTAGTGGATTTGGTGGAAGAACCATCTACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTAGTGCGC
210:ANTAGGGTGATACAGGACATGTATGACAGAGTCTCGACTAGCATTCAAACACCGGTAGGTTTGACAGAGACTTTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTAGTATGCCATCTTCTTCTTTGAAAAAAGG
221:CNAAACGAAGTTTTACTTTTAAGTATTCGGGTTAAAATTTATGGGCCTCGTGAATAGTATAGCCCAAAGTAATTGTGGGTTAAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTTACCATCTAGTACGCCGTCTTCTTCTTT

```
In general: i5 barcodes are in R2 reads; i7 barcodes in R1 reads. Sometimes, both of them are in the read *because the read might be smaller than 150 bp so both barcodes are being read*.


The barcode is not at the same place within the read - is that normal and expected? Sometimes 92, 114, 119 - 

Summary: 
I want to demultiplex paired-end reads using barcodes embedded within the read sequence, not the header. The barcode file uses i5/i7 columns per sample.


The reverse complement of the 5'-3' end PRIMER P5 (sequence: 'GTGTAGATCTCGGTGGTCGCCGTATCATTA' ) is present towards the end of R2 reads.
I could not find the normal sequence, reverse complement or complement sequence for PRIMER P7 (sequence: 'CAAGCAGAAGACGGCATACGAGAT') in my reads, R1 or R2. 



I got demuxbyname to work when using the *reverse-complement* of the i5 sequence. There is 50Gb that is unmatched still that I will try to assign using a higher Hamming distance. 










