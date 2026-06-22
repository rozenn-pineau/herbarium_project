
```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/jsmontgomery/anaconda/fastp

cd /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/demuxed_rc_twice/sequenced/

for r1 in *_R1.fastq.gz; do
        prefx=${r1%_R1.fastq.gz}

        fastp \
            --in1 ${prefx}_R1.fastq.gz \
            --in2 ${prefx}_R2.fastq.gz \
            --out1 ${prefx}_R1.unmerged.fq.gz \
            --out2 ${prefx}_R2.unmerged.fq.gz \
            --merge \
            --merged_out ${prefx}.collapsed.gz
done

```
