# BOUDREAU-PINSONNEAULT et al. SCENIC scRNA-seq analysis of reprogrammed Muller glia


We hope you find the explanations of the code easy to understand, we have included the explanations for a novice in bioinformatics as we once were when we started doing this analysis. 
We would love feedback on how we can improve and any issues you might be facing. Thank you for your time!

1) scRNA-seq RNA-velocity

We first used cellranger (the program created by 10x Genomics to align and count raw reads per cell with a genome) to generate the matrix files from our fastq files. For installation, follow this guide (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation). This will generate a folder containing the .bam files, matrix.mtx, barcodes.tsv and features.tsv files among others. 

Once you have this folder, the next step is to run Velocyto on your cellranger processed folder. We used the python version: http://velocyto.org/velocyto.py/tutorial/index.html. 

Side note: The hardest thing to setup are the installations/dependencies and 'environment' for your software (which in this case will be Velocyto). Since velocyto code is run on linux based environments, you would first have to setup linux on your pc. One way is to install ubuntu on your windows 10 pc using the official windows store: https://www.microsoft.com/en-ca/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab. Then you would have to install miniconda to set up the installation of velocyto. This is a blog outlining on how you could do that: https://varhowto.com/install-miniconda-ubuntu-20-04/. Once your conda environment is created, you can then run this code from the velocyto tutorial to start the installation of the dependencies:

conda install numpy scipy cython numba matplotlib scikit-learn h5py click

pip install velocyto

If you are not used to linux code, here are some basics:

if you want to look at your directory: 
ls

If you want to change directory: 
cd <new directory>

If you want to remove file:
rm <yourfiles.file>

If you want to remove folder:
rm -r <yourfolder>

We only had to use these commands for most of the code in linux. 

Since each machine is different, it would not be possible to be prepared for every possible error. That is why in case you run into an error, google is your friend! It will time and patience to find the solution to your specific problem but there will always be a solution on the web because your error won't be unique! 

Once installed, you can use this code to run velocyto on your 10x folder generated from cellranger:

velocyto run10x -m repeat_msk.gtf cellrangerfolder genes.gtf

Here, repeat_mask.gtf is the repeat mask files you can get from:

genes.gtf files you can get by running this code for the mouse genome:
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
Then unzip this file by using:
gunzip refdata-gex-mm10-2020-A.tar.gz
Then go to refdata-gex-mm10-2020-A/genes/genes.gtf for this file.

When velocyto is done, it will place a .loom file in your cellranger directory. This file will contain all the unspliced and spliced data from your raw reads per cells. We like working with this because it is smaller in file size.

The next step is the analysis itself. Now most of the work will be on python. To do this, we used Anaconda and downloaded Spyder through the Anaconda Navigator (https://docs.anaconda.com/anaconda/install/windows/). Amazing program that helps you use python and great for beginners. Once in spyder, use the python script called "Matrix plot final.py". This code will generate the UMAP and Matrix plots in the paper. 


2) SCENIC analysis pipeline

Luckily, if you have generated the .loom files from above, the next half is going to be fairly simple. The only thing I would recommend is having a powerful PC to run this next set of code because SCENIC is a very computationally intensive program. 

Firstly, install pyscenic on your spyder and linux environment because you need both to successfully generate the AUC UMAPs (https://pyscenic.readthedocs.io/en/latest/installation.html). 
  
Second, you need the auxillary datasets. These are detailed on the SCENIC page but for the analysis of our dataset, we used these three:
1)  https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl (These are motif annotations)
2)  https://github.com/aertslab/pySCENIC/blob/master/resources/mm_mgi_tfs.txt (These are names of the transcription factors)
3)  https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather (This file has TF network database)
  
Once you have these datasets and have pyscenic installed, use the python code called "SCENICIkzf14final.py" to run the initial steps for the SCENIC analysis. 
  
There are three lines that have to be run on a linux environment:
  
line  213: ./arboreto_with_multiprocessing.py mP0Scenic.loom mm_mgi_tfs.txt --method grnboost2 --output adjmP0.tsv --num_workers 10 --seed 777
  (Download arboreto_with_multiprocessing.py files from:  https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/cli/arboreto_with_multiprocessing.py)
  
line  227: pyscenic ctx adjtw.tsv mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname pbmc10k_filtered_scenic.loom --output reg.csv --mask_dropouts --num_workers 20
  
line  238: pyscenic aucell pbmc10k_filtered_scenic.loom reg.csv --output pyscenic_output.loom --num_workers 20
  
These lines will give you three important output files: adjmP0.tsv (adjacency matrix), reg.csv (regulon matrix), pyscenic_output.loom (integrated pyscenic loom file).
  
At the end of the code from "SCENICIkzf14final.py", you should be able to get an output loom file that has the AUC and louvain UMAPs embedded into it. 
  
For the next phase, use "SCENIC postanalysis.py". This will allow you to generate all the AUC regulons present in the manuscript. The final half of the code will give you the top5 regulon in each cluster.
  
This concludes the tutorial on how to perform SCENIC analysis on scRNA-seq datasets. Feel free to contact us for any questions. Thank you for reading, we hope this was helpful.
 
  
