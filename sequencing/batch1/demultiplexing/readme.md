This readme describes the steps from downloading the data from Novogene, demultiplexing and initial filtering steps/

## Retrieving the data from Novogene
### Downloading the data from Novogene using lftp
```
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC26058003-Z01-F002:uD2dCtY9@usftp22.novogene.com; mirror --verbose --use-pget-n=8 -c'
```
### Checking the download
```
md5sum -c MD5.txt
```

Every file looked good.

## Demultiplexing

(a) Match barcode to sequence

I used a custom R script to match the barcodes i5 and i7 to their corresponding sequence (match_barcode_sequence.R).

(b) demultiplex

demuxbyname from BBmap tools:


Using demuxbyname_batch1_barcodes.txt having i7 sequence then i5 sequence *reverse complemented* in the following format: NNNNNNNN+NNNNNNNN, Hamming distance of 1:


(This script took !3h with 367 Gbytes of data).

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /home/rozennpineau/java

R1=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/P1_WKDL260007250-1A_23C5CJLT4_L4_1.fq.gz
R2=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/P1_WKDL260007250-1A_23C5CJLT4_L4_2.fq.gz
barcodes=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxbyname_batch1_barcodes.txt
out_folder=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc

#move to bbmap tool folder
cd /home/rozennpineau/bbmap

./demuxbyname.sh in=$R1 in2=$R2 \
  out=$out_folder/sample_%_R1.fastq.gz out2=$out_folder/sample_%_R2.fastq.gz \
  delimiter=: prefixmode=f \
  names=$barcodes hdist=1 \
  outu=$out_folder/unmatched_R1.fastq.gz outu2=$out_folder/unmatched_R2.fastq.gz \
  -Xmx16g

```

It performed pretty well. 

The unmatch pile has 23G and 24Gbites:
```
23G Jun 19 00:42 unmatched_R2.fastq.gz 
24G Jun 19 00:42 unmatched_R1.fastq.gz
```
This gives file names that are barcode-based. I used a bash script to rename the files from their barcode combination to their sample names. 

script: rename_demux.sh


```
bash demuxed_rc/rename_demux.sh demuxbyname_batch1_barcode_to_sample.txt demuxed_rc/

==================================================
 Done.
   Renamed  : 508 files
   Missing  : 0 files (barcode not found in demuxed_rc/)
   Skipped  : 0 files (destination already existed)
==================================================

```


To get an idea of the distribution of the raw size between samples, I ran another script based on the sample size in Gbytes. 

script: get_raw_size_estimate.sh

```
bash demuxed_rc/get_raw_size_estimate.sh demuxed_rc/ sequenced_sample_list.txt sizes_summary.txt


==================================================
 Output written to : sizes_summary.txt
 Samples found     : 125
 Samples NOT found : 0  (listed as NOT FOUND in output)
==================================================


#All samples (not only those that were supposed to be sequenced
bash demuxed_rc/get_raw_size_estimate.sh demuxed_rc/ output_barcode_sequences.txt sizes_summary_all.txt

==================================================
 Output written to : sizes_summary_all.txt
 Samples found     : 254
 Samples NOT found : 1  (listed as NOT FOUND in output)
==================================================
```

I calculated the Hamming distance between each pair of barcodes that we used (distance within i5 barcodes, within i7 barcodes, then add the two). The minimum distance is 3 base pairs. So I am now running the demultiplexing script on the unmatched pile to see if we can recover some files from that data.

After running, the unmatch files are a bit smaller:
```
-rw-rw-r-- 1 rozennpineau rozennpineau   20G Jun 22 15:32 unmatched_R1.fastq.gz
-rw-rw-r-- 1 rozennpineau rozennpineau   19G Jun 22 15:32 unmatched_R2.fastq.gz

```

Renaming the files from barcode to file names:
```
bash rename_demux.sh demuxbyname_batch1_barcode_to_sample.txt demuxed_rc/unmatched_pile/

```

Concatenate the files with the same names in a third folder: 
```
for f in /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc/unmatched_pile/*.fastq.gz; do
    fname=$(basename "$f")
    cat "$f" "/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc/$fname" > "demuxed_rc_twice/$fname"
done

# Check in every folder that the sun of the two equals the new file using
zcat sample_164_R1.fastq.gz | wc -l
# worked
```
Barcodes correspond to samples that were not sequences (but prepped). I want now to focus only on the samples that were sequenced and I want to look at how big the samples that are not supposed to be there are.

```
mkdir -p sequenced not_sequenced

# Move files matching the sample list
for id in $(cat ../sequenced_sample_list.txt); do
    mv sample_${id}_R1.fastq.gz sequenced/ 2>/dev/null
    mv sample_${id}_R2.fastq.gz sequenced/ 2>/dev/null
done

# Move everything remaining
mv sample_*.fastq.gz not_sequenced/ 2>/dev/null

``` 

25 GBytes of samples that were not supposed to be sequenced...

304 GBytes of samples that were supposed to be sequenced.



