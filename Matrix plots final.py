# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:01:52 2021

@author: javeda
"""
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

new_cluster_names = [
    'Control Müller-1', 'Control Cones',
    'Control Rods','Control Müller-2', 'Control Bipolar cells', 'Control Endothelial cells', 'Control Microglia']
adata.rename_categories('leiden', new_cluster_names)
#You can find genes enriched per cluster with 3 different statistical tests.
Control = adata



#Change directory
wdir = "D:/10x/Casz1/Muller/Ikzf14"
os.chdir( wdir )
#Your loom file is read here. If you have set the working directory as the folder that contains your loom file, then just add the name of the .loom file here.
adata = scv.read('IK14.loom', cache=True)
adata.var_names_make_unique()

#Parameters to sort cells based on genes expressed, play around to see what your data looks like.
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#Parameters to sort cells based on ncount and mito counts, another parameter to play around with.
adata
adata.var['M'] = adata.var_names.str.startswith('M')  # annotate the group of mitochondrial genes as 'm'
sc.pp.calculate_qc_metrics(adata, qc_vars=['M'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_M'],
             jitter=0.4, multi_panel=True)


#normalization using scVelo
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=5000)

#Gene enriched in the cells and PCA/UMAP clustering of the data.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=5, n_pcs=40)
sc.tl.umap(adata, min_dist=(0.5))
#Leiden clusters your data. Everything referred to as 'leiden' in the rest of the code are your clusters.
sc.tl.leiden(adata, resolution=0.15)
sc.pl.umap(adata, color=['leiden'])

new_cluster_names = [
    'Ikzf1/4 Müller-1', 'Ikzf1/4 Rods',
    'Ikzf1/4 Bipolar cells 1','Ikzf1/4 Bipolar cells 2', 'Ikzf1/4 Microglia', 'Ikzf1/4 Cones', 'Ikzf1/4 Bipolar cells 3', 'Ikzf1/4 Müller-2', 'Ikzf1/4 Bipolar cells 4', 'Ikzf1/4 Bipolar cells 5', 'Ikzf1/4 Endothelial cells-1', 'Ikzf1/4 Endothelial cells-2', 'Ikzf1/4 Pericytes', 'Ikzf1/4 Müller-3']
adata.rename_categories('leiden', new_cluster_names)
#You can find genes enriched per cluster with 3 different statistical tests.
Ikzf14 = adata
    


#Change directory
wdir = "D:/10x/Casz1/mP14"
os.chdir( wdir )
#Your loom file is read here. If you have set the working directory as the folder that contains your loom file, then just add the name of the .loom file here.
adata = scv.read('P14_CR.loom', cache=True)
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
scv.pp.filter_and_normalize(adata, min_shared_counts=3, n_top_genes=25000)

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
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])


new_cluster_names = [
    '0', 'P14 Bipolar cells',
    '2','3', '4', '5', '6', '7']
adata.rename_categories('leiden', new_cluster_names)

P14 = adata

outer = Ikzf14.concatenate([Control], join='outer')
outer

    
marker_genes_dictBP = {'BC1A': ['Pcdh17'],'BC1B': ['Pcdh10'],
    'BC2/BC5B': ['Chrm2'],
    'BC3A': ['Erbb4'],'BC3B': ['Nnat'],
    'BC4': [ 'Col11a1'],
    'BC5A': ['Sox6'], 'BC5D': ['Lrrtm1'], 'BC6': ['Cck'], 'BC7': ['Igfn1'],'BC8/9': ['Serpini1','Cpne9'], 'RBC': ['Vstm2b','Casp7'], 'Cones' : ['Thrb']}

marker_genes_dict = {'Pr markers': ['Crx', 'Rcvrn', 'Pdc'],
    'Rods': [ 'Rho','Pde6b'],
    'Cones': ['Arr3', 'Opn1sw','Opn1mw'],
    'Bipolars': [ 'Vsx2', 'Grm6','Otx2'],
    'Müller glia': ['Clu', 'Aqp4', 'Rlbp1', 'Slc1a3']}

marker_genes_dictyap = {
    'Bipolar markers': ['Isl1', 'Cabp5', 'Grm6', 'Neurod4','Pcp2']}

BP = outer[outer.obs['leiden'].isin(['Control Bipolar cells',
    'Ikzf1/4 Bip/Co 1','Ikzf1/4 Bip/Co 2', 'Ikzf1/4 Bip/Co 3', 'Ikzf1/4 Bip/Co 4', 'Ikzf1/4 Bip/Co 5'])]
All = outer[outer.obs['leiden'].isin([ 'Control Müller-1', 'Control Cones',
    'Control Rods','Control Müller-2', 'Control Bipolar cells', 'Ikzf1/4 Müller-1', 'Ikzf1/4 Rods',
    'Ikzf1/4 Bip/Co 1','Ikzf1/4 Bip/Co 2',  'Ikzf1/4 Cones', 'Ikzf1/4 Bip/Co 3', 'Ikzf1/4 Müller-2', 'Ikzf1/4 Bip/Co 4', 'Ikzf1/4 Bip/Co 5',  'Ikzf1/4 Müller-3'])]
mBP = outer[outer.obs['leiden'].isin([    'Ikzf1/4 Bip/Co 1','Ikzf1/4 Bip/Co 2', 'Ikzf1/4 Bip/Co 3', 'Ikzf1/4 Bip/Co 4', 'Ikzf1/4 Bip/Co 5'])]

adata.var_names = adata.var_names

ax = sc.pl.heatmap(BP, marker_genes_dictBP, groupby='leiden', cmap='viridis', dendrogram=False)
sc.pl.dotplot(outer, marker_genes_dictyap, groupby='leiden', dendrogram=False)
sc.pl.matrixplot(mBP, marker_genes_dictBP, 'leiden', dendrogram=False, standard_scale='var')
sc.pl.matrixplot(outer, marker_genes_dictyap,'leiden', dendrogram=False, standard_scale='var')





