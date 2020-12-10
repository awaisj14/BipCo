# BOUDREAU-PINSONNEAULT et al. SCENIC and scVelo scRNA-seq analysis of reprogrammed Muller glia


We hope you find the explanations of the code easy to understand, we have included the explanations for a novice in bioinformatics as we were once when we started doing this analysis. 
We would love feedback on how we can improve and any issues you might be facing. Thank you for your time!

1) scRNA-seq RNA-velocity

We first used cellranger (the program created by 10x Genomics to align and count raw reads per cell with a genome) to generate the matrix files from our fastq files. This will generate a folder containing the .bam files, matrix.mtx, barcodes.tsv and features.tsv files among others. 

Once you have this folder, the next step is to run RNA-velocyto on your cellranger processed folder. We used the python version: http://velocyto.org/velocyto.py/tutorial/index.html. 
