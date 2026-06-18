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











