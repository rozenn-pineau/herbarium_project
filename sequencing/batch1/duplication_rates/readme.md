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
        output=$out/${name}.prefixed.bam

        #check if file is present
        if [ -s "$output" ]; then
                echo "$output already exists and is not empty, skipping."
        else
        samtools view -h ${name}.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/M_/' | samtools view -bS - > $out/${name}.prefixed.bam
        fi

        #sort prefixed bam
        sorted_output=$out/${name}.prefixed.sorted.bam
        if [ -s "$sorted_output" ] ; then
                echo "$sorted_output already exists and is not empty, skipping."
        else
                #sort and index
                sambamba sort -m 15GB --tmpdir tmp -t $threads -o $sorted_output/${name} $output
        fi
done

#checking for the presence of the file (as compared to the completion) works in this scenario because samtools outputs only when it is done, not as it goes
```
### Adding F_ and R_ to forward and reverse reads

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba
module load samtools

### Adding F_ to forward reads

cd /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged
mkdir -p forward #creates the folder if it is missing, ignores it safely if it exists
threads=6

for bam in *.unmerged.sorted.bam; do

        name=${bam%.unmerged.sorted.bam}

        if [ -s "forward/${name}.F.prefixed.unmerged.bam" ] ; then
                echo "File forward/${name}.F.prefixed.unmerged.bam exists, skipping."
        else
                # split
                samtools view -F 0x10 -h ${name}.unmerged.sorted.bam | samtools view -bS - > forward/${name}.F.unmerged.sorted.bam #exclude reverse read
                # rename
                samtools view -h forward/${name}.F.unmerged.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/F_/' | samtools view -bS - > forward/${name}.F.prefixed.unmerged.bam
        fi

        #sort and index
        echo "Sorting forward/${name}.F.prefixed.unmerged.bam."
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o forward/${name}.F.prefixed.unmerged.sorted.bam forward/${name}.F.prefixed.unmerged.bam

done


### Adding R_ to reverse reads

cd /scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged
mkdir -p reverse
threads=6

for bam in *.unmerged.sorted.bam; do

        name=${bam%.unmerged.sorted.bam}

        if [ -s "reverse/${name}.R.prefixed.unmerged.bam" ]; then
                echo "File reverse/${name}.R.prefixed.unmerged.bam exists, skipping."
        else
                # split
                samtools view -f 0x10 -h ${name}.unmerged.sorted.bam | samtools view -bS - > reverse/${name}.R.unmerged.sorted.bam #exclude forward read
                # rename
                samtools view -h reverse/${name}.R.unmerged.sorted.bam | sed 's/^@/&/;/^[^@]/s/^/R_/' | samtools view -bS - > reverse/${name}.R.prefixed.unmerged.bam
        fi

        #sort and index
        echo "Sorting reverse/${name}.R.prefixed.unmerged.bam."
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o reverse/${name}.R.prefixed.unmerged.sorted.bam reverse/${name}.R.prefixed.unmerged.bam

done

```


### Merge prefixed bams

!!! Add a sorting step after the merging step !

```
module load samtools
unmerged_R_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged/reverse
unmerged_F_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged/forward
collapsed_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/collapsed/renamed
output_folder=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/final

cd $unmerged_R_bams

for r1 in *.R.prefixed.unmerged.sorted.bam; do

        prefx=${r1%.R.prefixed.unmerged.sorted.bam}
        samtools merge $output_folder/${prefx}.prefixed.bam $unmerged_R_bams/${prefx}.R.prefixed.unmerged.sorted.bam $unmerged_F_bams/${prefx}.F.prefixed.unmerged.sorted.bam $collapsed_bams/${prefx}.prefixed.sorted.bam

done


```

### Coverage
Estimate coverage by calculting the number of mapped bases per sample.

```
module load samtools
path=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/final
cd $path

out=number_base_pairs_per_sample_.txt
echo -e "sample\tnum_base_pairs" > $out

for bam in *.prefixed.bam; do

  name=${bam%.prefixed.bam}
  echo "Calculating base pairs mapped in $name..."
  num=$(samtools stats $bam | grep "bases mapped (cigar):" | cut -f 3)
  echo -e "$name\t$num" >> $out

done

```

### Split bams into scaffolds

I first have to split the bam into chromosomes, as it keeps running out of memory on the full bam file.

#create a subfolder for each file
#split bam into chromosomes

```

cd /scratch/midway2/rozennpineau/herbarium/bams/final/sorted

module load samtools

for bam in *.prefixed.sorted.bam; do
        
        name=${bam%.prefixed.sorted.bam}
       
        
        #check if sample has already been split
        if [ -s "split/$name/scaffold_999.bam.bai" ]; then
                echo "$name has been split, moving on."
        else
                 mkdir -p split/$name
                #retrieve scaffolds with samtools idxstats, ignore last line of output as it is a wildcard
                for chr in $(samtools idxstats $bam | cut -f1 | head -n -1); do

                        echo "Splitting $chr..."
                        samtools view -b $bam "$chr" > split/$name/"${chr}.bam" #split into scaffolds
                        samtools index split/$name/"${chr}.bam" #index

                done
        fi

done

```

### Sample 29 was left behind - troubleshooting when it failed

Troubleshooting sample 29
Reverse reads: 
samtools view sample_29.R.unmerged.sorted.bam | tail - OK
Forward reads
samtools view sample_29.F.unmerged.sorted.bam | tail - OK
Collapsed reads
samtools view sample_29.prefixed.sorted.bam | tail
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
— not OK

The step before sorting: 
samtools view sample_29.prefixed.bam | tail -- looks good ! re-sort for here

```
#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/sambamba

module load samtools
unmerged_R_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged/reverse
unmerged_F_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/unmerged/forward
collapsed_bams=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/sorted/collapsed/renamed
output_folder=/scratch/midway3/rozennpineau/herbarium/01.RawData/P1/bams/final

threads=6

#sort sample 29
cd $collapsed_bams
sambamba sort -m 15GB --tmpdir tmp -t $threads -o sample_29.prefixed.sorted.bam sample_29.prefixed.bam

cd $unmerged_R_bams

for r1 in sample_29.R.prefixed.unmerged.sorted.bam; do

        prefx=${r1%.R.prefixed.unmerged.sorted.bam}
        samtools merge -f $output_folder/${prefx}.prefixed.bam $unmerged_R_bams/${prefx}.R.prefixed.unmerged.sorted.bam $unmerged_F_bams/${prefx}.F.prefixed.unmerged.sorted.bam $collapsed_bams/${prefx}.prefixed.sorted.bam
        #-f : allow overwriting

        #sort
        sambamba sort -m 15GB --tmpdir tmp -t $threads -o $output_folder/sorted/${prefx}.prefixed.sorted.bam $output_folder/${prefx}.prefixed.bam

done
```

These all worked, going back to "splitting the bam into scaffolds step". 

### Running Dedup


!!! Add a way to follow how many files have been processed because it is hard to read the slurm.err or .out and know the progression of the work !!!

```
#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --output=dedup.out
#SBATCH --error=dedup.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20GB

module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/dedup/

module load samtools

working_dir=/scratch/midway2/rozennpineau/herbarium/01.RawData/P1/bams/final/sorted/split
cd $working_dir

for dir in ./sample_*; do
    cd $dir
        for bam in *.bam; do
                name=${bam%.bam}
                mkdir -p $name
                dedup -i $bam -o ${name}
        done
    cd ../
done
```

This script was running very slow, so I updated it with some parallelization: 

I also realized I made a mistake and the bams were being overwritten and the log was not properly stored, so I deleted the splitted bams to run that script again. 

```
#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --output=dedup.out
#SBATCH --error=dedup.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=8GB


module load python/anaconda-2021.05
source /software/python-anaconda-2021.05-el7-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/dedup/
module load samtools
module load parallel   # load GNU parallel if it's provided as a module on Midway2

working_dir=/scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/to_do
cd "$working_dir"
output_folder=/scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup

#find: lists all bams in the directories
find . -mindepth 2 -maxdepth 2 -name "*.bam" | sed 's|^\./||' | \
#-mindepth 2 = only look at files that are 2 levels deep 
#-maxdepth 2 = don't go any deeper than that 

#running dedup on each file, in parallel
# -j N = number of dedup jobs to run at once; match to --cpus-per-task above
parallel -j 8 --joblog parallel_dedup.log '
    dir=$(dirname {})
    fname=$(basename {})
    name=${fname%.bam}
    cd "$dir"
    mkdir -p "$output_folder/$name"
    dedup -i "$fname" -o "$output_folder/$name"
'
```

### Retrieve duplication rate from files after running DeDup

```
#!/bin/bash
#SBATCH --job-name=calc_dup_rate
#SBATCH --output=calc_dup_rate.out
#SBATCH --error=calc_dup_rate.err
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12GB

cd /scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup

echo -e "sample\tmean duplication rate" > mean_dup_rate.txt
for sample_dir in sample*/; do
cd $sample_dir
  echo -e "scaffold\tduplication rate" > dup_rate_summary.txt
  for scaffold_dir in [Ss]*/; do
    scaffold_name=$(basename "$scaffold_dir")
    
    dup_rate=$(awk '/Duplication Rate/ {print $3}' "$scaffold_dir/$scaffold_name.log") #Retrieve duplication rate

    echo -e "${scaffold_name}\t${dup_rate}" >> dup_rate_summary.txt
  done
  #calculate mean over all scaffolds
  mean=$(awk 'NR>1 {sum += $2; n++} END {print sum/n}' dup_rate_summary.txt)
cd ..
echo -e "$sample_dir\t$mean" >> mean_dup_rate.txt
done
```

Calculate mean duplication rate for Scaffolds 1-16:

```
echo -e "Sample\tScaffold1-16_duplication_rate" > scaffold1-16_mean_dup_rate.txt

for dir in ./sample*; do
  #echo $dir
  dup_rate=$(cat $dir/dup_rate_summary.txt | grep Scaffold | awk 'NR>1 {sum += $2; n++} END {print sum/n}')
  echo -e "${dir}\t${dup_rate}" >> Scaffold1-16_mean_dup_rate.txt
done
```

### Merge bams after Dedup

```
start_dir=/scratch/midway2/rozennpineau/herbarium/bams/final/sorted/split/dedup

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bamtools
module load samtools

ulimit -n 4096 # increase upper limit of number of files that can be opened at once

cd $start_dir
for dir in ./*; do #list directories one level down only 
    cd $dir

    realpath */[Ss]*rmdup.bam > bams_to_merge.list
    echo "merging $dir..."
    bamtools merge -list bams_to_merge.list -out $dir.scaffolds.dedup.bam
    samtools index $dir.scaffolds.dedup.bam #index

    cd ..
done
```