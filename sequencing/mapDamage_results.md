### Understanding what mapDamage scores and what to look at. 

link to publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC3694634/pdf/btt193.pdf


How to understand the key parameters being estimated and the key files being generated. 

Parameter: Description
λ (Lambda):	The probability of terminating an overhang.

D (Delta D):	Cytosine deamination probability in a double-stranded context.

S (Delta S):	Cytosine deamination probability in a single-stranded context.

θ (Theta):	The mean difference rate between the reference and sample that is not due to DNA damage (e.g., natural polymorphism).

v (Nu):	Nick frequency within the DNA fragment.


Summary of Key Output Files
File Name: Content Description
3pGtoA_freq.txt:	Frequencies of G → A substitutions for the first N bases from the 3' end.

5pCtoT_freq.txt:	Frequencies of C → T substitutions for the first N bases from the 5' end.

dnacomp.txt:	Counts of each base (A, T, G, C) across all reads and surrounding reference regions.

lgdistribution.txt:	Read length distributions according to the strand (+/-) they map to.

alignment.rescaled.bam:	Generated if the --rescale option is used. It downscales quality scores of bases likely affected by deamination to reduce noise in downstream analysis like SNP calling.
