### Understanding the output from Picard Mark Duplicates

Picard Mark Duplicates outputs a txt file in addition to the bam marked duplicates. 

The txt file has several entries:

UNPAIRED_READS_EXAMINED
--> Single-end reads (or orphan reads) that were checked for duplicates.


READ_PAIRS_EXAMINED
--> Properly paired reads examined (this is the bulk of your data)


SECONDARY_OR_SUPPLEMENTARY_RDS
→ Reads that are:
secondary alignments (multi-mapping)
supplementary alignments (split reads)

--> These are ignored for duplicate marking, but useful to know if you have lots of multi-mapping (common in repetitive genomes)


UNMAPPED_READS
--> Reads that didn’t align to your reference


UNPAIRED_READ_DUPLICATES
--> Duplicate reads among the unpaired reads


READ_PAIR_DUPLICATES
--> Duplicate paired-end fragments


READ_PAIR_OPTICAL_DUPLICATES
--> Subset of duplicates likely due to optical artifacts (same cluster detected multiple times on the sequencer)


ESTIMATED_LIBRARY_SIZE
--> This estimates how many unique DNA fragments were in your original library.

I will look at those different proportions and how they relate to our library preparation metrics using R.
















