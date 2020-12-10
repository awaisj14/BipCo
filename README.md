# BOUDREAU-PINSONNEAULT et al. SCENIC and scVelo scRNA-seq analysis of reprogrammed Muller glia


We hope you find the explanations of the code easy to understand, we have included the explanations for a novice in bioinformatics as we were once when we started doing this analysis. 
We would love feedback on how we can improve and any issues you might be facing. Thank you for your time!

1) scRNA-seq RNA-velocity

We first used cellranger (the program created by 10x Genomics to align and count raw reads per cell with a genome) to generate the matrix files from our fastq files. This will generate a folder containing the .bam files, matrix.mtx, barcodes.tsv and features.tsv files among others. 

Once you have this folder, the next step is to run Velocyto on your cellranger processed folder. We used the python version: http://velocyto.org/velocyto.py/tutorial/index.html. 

Side note: The hardest thing to setup are the installations/dependencies and 'environment' for your software (which in this case will be Velocyto). Since velocyto is run on linux based environments, you would first have to setup linux on your pc. One way is to install ubuntu on your windows 10 pc using the official windows store: https://www.microsoft.com/en-ca/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab. Then you would have to install miniconda to set up the installation of velocyto. This is a blog outlining on how you could do that: https://varhowto.com/install-miniconda-ubuntu-20-04/. Once your conda environment is created, you can then run this code from the velocyto tutorial to start the installation of the dependencies:

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

The next step is the analysis itself. Now most of the work will be on python. To do this, we used Spyder https://www.spyder-ide.org/ using Anaconda. Amazing program that helps you use python and great for beginners. Install it and then you are not in the python environment. We use scVelo for the next analysis step. Just like linux, setting up the dependencies is the first step. We stitched together code from different places to make this python code: https://github.com/awaisj14/BipCo/blob/main/Mullers%20scVelo.py

There you have it, you should be able to visualize your gene of interest in your dataset.

2) SCENIC analysis pipeline

Luckily, if you have generated the .loom files from above, the next half is going to be fairly simple. The only thing i would recommend is having a powerful PC to run this next set of code because SCENIC is very computationally intensive program. 

