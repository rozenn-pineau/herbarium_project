Sequencing batch 2 pre run: 540 samples pooled into 34 pools, sequenced on 2 lane using NovaSeq 25B. The goal is to have an idea of how much data we get when we sequence those samples together. Equimolarity does not equate equal sequence coverage out of the sequencer. We sequenced all 540 samples at low coverage (~1.7x per sample) to recalibrate our pool calculations for deeper sequencing. 

This readme describes the steps from downloading the data from Novogene, demultiplexing and initial filtering steps.

## Retrieving the data from Novogene
# Downloading the data from Novogene using lftp

```
conda install conda-forge::lftp
conda activate /project/kreiner/rpineau/lftp

lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC26087787-Z01-F001:OHmE9E4M@usftp22.novogene.com; mirror --verbose --use-pget-n=8 -c'
```

### Saving the data on the NAS
(1) Ssh to the NAS: ssh rpineau@kreinerlab.uchicago.edu (make sure to log onto Cisco first)
(2) navigate where you want to save the data (/volume1/Data_2025/herbarium/20260903)
(3) 
```
/volume1/Data_2025/herbarium/20260903
```

### Checking the download
```
P18_L2_WKDL260014409-1A_23C557LT4_L7_1.fq.gz: OK
P18_L2_WKDL260014409-1A_23C557LT4_L7_2.fq.gz: OK
P1_L1_WKDL260014392-1A_23C557LT4_L8_1.fq.gz: OK
P1_L1_WKDL260014392-1A_23C557LT4_L8_2.fq.gz: OK
```

Every file looked good.

## Demultiplexing

To note: 
I forgot that a few samples had the same barcode combination when optimizing the pooling strategy, and some of them ended up on the same lane (I won't be able to separate them), while some of them ended up on different lanes. For those on different lanes, I can recover the identity of those samples if I demultiplex the lanes separately. 

(a) Match barcode to sequence

I used a custom R script to match the barcodes i5 and i7 to their corresponding sequence (match_barcode_sequence_batch2.R, also in this folder).
On NovaSeq the read 2 read goes in the opposite direction, so the instrument actually sequences through i5 from the other end --> the i5 sequences recorded in the FASTQ read headers will be the reverse complement of what we designed.

!! Use the *reverse complement sequence for i5* !!



(b) demultiplex

demuxbyname from BBmap tools:


Using demuxbyname_batch1_barcodes.txt having i7 sequence then i5 sequence *reverse complemented* in the following format: NNNNNNNN+NNNNNNNN, Hamming distance of 1:


(This script took 4 hours, 45 minutes for lane 1, 6h 20 min for lane 2).

```
#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --output=demux_l1.out
#SBATCH --error=demux_l1.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=16G   # memory per cpu-core

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /home/rozennpineau/java

R1=/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/P18_L2_WKDL260014409-1A_23C557LT4_L7_1.fq.gz
R2=/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/P18_L2_WKDL260014409-1A_23C557LT4_L7_2.fq.gz
barcodes=/scratch/midway3/rozennpineau/herbarium/batch2/barcodes/lane2_demuxbyname_batch2_barcodes.txt
out_folder=/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/demuxed

mkdir -p $out_folder

#move to bbmap tool folder
cd /home/rozennpineau/bbmap

#! make sure to adjust the Hamming distance to the value that matters to you!

./demuxbyname.sh in=$R1 in2=$R2 \
  out=$out_folder/sample_%_R1.fastq.gz out2=$out_folder/sample_%_R2.fastq.gz \
  delimiter=: prefixmode=f \
  names=$barcodes hdist=1 \
  outu=$out_folder/unmatched_R1.fastq.gz outu2=$out_folder/unmatched_R2.fastq.gz \
  threads=4 \
  -Xmx60g

```

Some stats:


Lane 2 data (/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/demuxed): 

total 522G
  26G -rw-rw-r-- 1 rozennpineau rozennpineau   26G Sep  5 17:26 unmatched_R2.fastq.gz
  27G -rw-rw-r-- 1 rozennpineau rozennpineau   27G Sep  5 17:26 unmatched_R1.fastq.gz

Input is being processed as paired
Time:               22100.241 seconds.
Reads Processed:    8429576342  381.42k reads/sec
Bases Processed:    1264436451300       57.21m bases/sec
Reads Out:          7721170316
Bases Out:          1158175547400
Yield:              0.91596

Lane 1 data (/scratch/midway2/rozennpineau/herbarium/batch2/raw/P1_L1/demuxed):

total 470G
  42G -rw-rw-r-- 1 rozennpineau rozennpineau   42G Sep  5 16:16 unmatched_R2.fastq.gz
  44G -rw-rw-r-- 1 rozennpineau rozennpineau   44G Sep  5 16:16 unmatched_R1.fastq.gz
Time:               17109.696 seconds.
Reads Processed:    7446614190  435.23k reads/sec
Bases Processed:    1116992128500       65.28m bases/sec
Reads Out:          6249535264
Bases Out:          937430289600
Yield:              0.83925



This gives file names that are barcode-based. I used a bash script to rename the files from their barcode combination to their sample names.
(script is very fast and renames the files in place)

Lane 1 data is on Midway2 - 
```
barcodes=/scratch/midway3/rozennpineau/herbarium/batch2/barcodes/lane1_demuxbyname_batch2_barcode_to_sample.txt
demuxed_folder=/scratch/midway2/rozennpineau/herbarium/batch2/raw/P1_L1/demuxed
bash /scratch/midway3/rozennpineau/herbarium/scripts/rename_demux.sh $barcodes $demuxed_folder

==================================================
 Done.
   Renamed  : 540 files
   Missing  : 0 files (barcode not found in /scratch/midway2/rozennpineau/herbarium/batch2/raw/P1_L1/demuxed)
   Skipped  : 0 files (destination already existed)
==================================================
```

How is the data distributed between samples?
```
scripts_folder=/scratch/midway3/rozennpineau/herbarium/scripts
demuxed_folder=/scratch/midway2/rozennpineau/herbarium/batch2/raw/P1_L1/demuxed
sample_list=/scratch/midway3/rozennpineau/herbarium/batch2/barcodes/lane1_demuxbyname_batch2_barcode_to_sample.txt # list of samples in column 2

bash $scripts_folder/get_raw_size_estimate.sh $demuxed_folder $sample_list lane1_raw_size_summary.txt 2

==================================================
 Output written to : lane1_raw_size_summary.txt
 Samples found     : 270
 Samples NOT found : 0  (listed as NOT FOUND in output)
==================================================

```



**Check that the conversion is correct?**

Lane 2 data is on Midway3 - 
```
barcodes=/scratch/midway3/rozennpineau/herbarium/batch2/barcodes/lane2_demuxbyname_batch2_barcode_to_sample.txt
demuxed_folder=/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/demuxed
bash /scratch/midway3/rozennpineau/herbarium/scripts/rename_demux.sh $barcodes $demuxed_folder

 Done.
   Renamed  : 532 files
   Missing  : 8 files (barcode not found in /scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/demuxed)
   Skipped  : 0 files (destination already existed)


   WARN: not found -> sample_ACGTTACC+AAGTGTCG_R1.fastq.gz
WARN: not found -> sample_ACGTTACC+AAGTGTCG_R2.fastq.gz
OK  : sample_ACGTTACC+CACAAGTC_R1.fastq.gz  ->  sample_617_R1.fastq.gz
OK  : sample_ACGTTACC+CACAAGTC_R2.fastq.gz  ->  sample_617_R2.fastq.gz
WARN: not found -> sample_ACGTTACC+AGTCTCAC_R1.fastq.gz
WARN: not found -> sample_ACGTTACC+AGTCTCAC_R2.fastq.gz
WARN: not found -> sample_ACGTTACC+CATGGAAC_R1.fastq.gz
WARN: not found -> sample_ACGTTACC+CATGGAAC_R2.fastq.gz
WARN: not found -> sample_ACGTTACC+CTCAGCTA_R1.fastq.gz
WARN: not found -> sample_ACGTTACC+CTCAGCTA_R2.fastq.gz

```
How is the data distributed between samples?
```
scripts_folder=/scratch/midway3/rozennpineau/herbarium/scripts
demuxed_folder=/scratch/midway3/rozennpineau/herbarium/batch2/raw/01.RawData/P18_L2/demuxed
sample_list=/scratch/midway3/rozennpineau/herbarium/batch2/barcodes/lane2_demuxbyname_batch2_barcode_to_sample.txt # list of samples in column 2

bash $scripts_folder/get_raw_size_estimate.sh $demuxed_folder $sample_list lane2_raw_size_summary.txt 2

==================================================
 Output written to : lane2_raw_size_summary.txt
 Samples found     : 266
 Samples NOT found : 4  (listed as NOT FOUND in output)
==================================================

```

Unmatched stats

Lane 2
26G -rw-rw-r-- 1 rozennpineau rozennpineau 26G Sep  5 17:26 unmatched_R2.fastq.gz
27G -rw-rw-r-- 1 rozennpineau rozennpineau 27G Sep  5 17:26 unmatched_R1.fastq.gz


## Checking the quality with fastp

Parallelized version of the script.

```
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=10G

# activate conda
module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/jsmontgomery/anaconda/fastp

WORKDIR=/scratch/midway2/rozennpineau/herbarium/batch2/demuxed
cd $WORKDIR

run_fastp() {
    prefx=${1%_R1.fastq.gz}

    fastp \
        --in1 ${prefx}_R1.fastq.gz \
        --in2 ${prefx}_R2.fastq.gz \
        --out1 ${prefx}_R1.unmerged.fq.gz \
        --out2 ${prefx}_R2.unmerged.fq.gz \
        --merge \
        --merged_out ${prefx}.collapsed.fq.gz \
        --html ${prefx}.html \
        --json ${prefx}.json \
        --thread 2
}

export -f run_fastp

parallel -j $SLURM_NTASKS_PER_NODE run_fastp ::: *_R1.fastq.gz


```