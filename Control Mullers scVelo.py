# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:33:05 2020

@author: javeda
"""
#Installation of libraries and their dependencies
pip install -U scvelo
pip install anndata scanpy numpy scipy pandas scikit-learn matplotlib
pip install python-igraph louvain
pip install python-igraph --upgrade --quiet


#Import function is similar to the library() function in R.
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import loompy as lp
import matplotlib as mpl

#Settings to beautify the images
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=600, facecolor='white')

#Change directory
wdir = "D:/10x/Casz1/Muller/Control"
os.chdir( wdir )
#Your loom file is read here. If you have set the working directory as the folder that contains your loom file, then just add the name of the .loom file here.
adata = scv.read('Control.loom', cache=True)
adata.var_names_make_unique()

#Parameters to sort cells based on genes expressed, play around to see what your data looks like.
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#Parameters to sort cells based on ncount and mito counts, another parameter to play around with.
adata
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


#normalization using scVelo
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)

#Gene enriched in the cells and PCA/UMAP clustering of the data.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=40)
sc.tl.umap(adata, min_dist=(0.5))
#Leiden clusters your data. Everything referred to as 'leiden' in the rest of the code are your clusters.
sc.tl.leiden(adata, resolution=0.125)
sc.pl.umap(adata, color=['leiden'])

sc.pl.umap(adata, color='leiden', legend_loc='on data',legend_fontsize='8', title='Control', frameon=False, size=15, save='Controlondata.pdf')
    Control = adata
new_cluster_names = [
    'Müller-1', 'Cones',
    'Rods','Müller-2', 'Bipolar cells', 'Endothelial cells', 'Microglia']
adata.rename_categories('leiden', new_cluster_names)
#You can find genes enriched per cluster with 3 different statistical tests.

sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#scVelo - Velocity calculation time!
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', save='Controltrajectory.pdf')
scv.pl.velocity_embedding(adata, arrow_length=5, arrow_size=1, dpi=600)

#Here you can look at the spliced and unspliced concentration of your genes. Code for multiple.
scv.pl.velocity(adata, ['Ikzf1', 'Thrb', 'Pcp2'], ncols=3)

#....and code for single genes.
scv.pl.velocity(adata, 'Null', save='Ikzf1velo.pdf')
scv.pl.velocity(adata, 'Thrb', save='Thrbvelo.pdf')
#Using scVelo for ranking your genes based on clusters.
scv.tl.rank_velocity_genes(adata, groupby='leiden', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

#Calculating Velocity length and confidence here.
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

#Cycling progenitors
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

#Calculating velocity trajectories here.
scv.pl.velocity_graph(adata, threshold=.1)
