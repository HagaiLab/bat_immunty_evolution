import numpy as np
import scanpy as sc
import os
import pandas as pd
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt

"""
This script clusters and plotts the expression of top 100 genes from each state among all 3-states cells.
Plotting is done using sns.clustermap 
"""

# Load h5ad data
bat1k_all_data = sc.read(r"C:\Users\TzachiHNB6\Documents\annotated_count_matrices\bat\blood\bat_mouse_integration\PBMCs_annotated_raw_bat1k\PBMCs_annotated_raw_bat1k.h5ad")
bat_metadata = pd.read_csv(r"C:\Users\TzachiHNB6\Downloads\pbmcs_bat1k_metadata.csv")
bat1k_all_data.obs = bat_metadata
bat1k_all_data.obs.rename(columns={"Unnamed: 0": "cell"}, inplace=True)
bat1k_all_data.obs.set_index("cell", inplace=True)

# Keep only CD16 + IFN+ Monocytes in Unstimulated and Poly(I:C)
bat1k_all_data_cd16_ifn_pic = bat1k_all_data[(bat1k_all_data.obs.batch.isin(['Poly(I:C)',"Lipofectamine",'Unstimulated'])) & (bat1k_all_data.obs.inter_species_ann.isin(['Monocytes_CD16', "Monocytes_IFNB1"]))&(bat1k_all_data.obs.Sample!='Bat2')&(bat1k_all_data.obs.pct_counts_MT<10)]

# Normalize and scale data
sc.pp.normalize_total(bat1k_all_data_cd16_ifn_pic)
sc.pp.log1p(bat1k_all_data_cd16_ifn_pic)
sc.pp.scale(bat1k_all_data_cd16_ifn_pic, max_value=10)


# Add states to data
conditions = [
    (bat1k_all_data_cd16_ifn_pic.obs['inter_species_ann'] == 'Monocytes_IFNB1'),
    (bat1k_all_data_cd16_ifn_pic.obs['inter_species_ann'] == 'Monocytes_CD16') & 
    (bat1k_all_data_cd16_ifn_pic.obs['batch'] == 'Poly(I:C)')]
choices = ['IFN+', 'IFN-']
bat1k_all_data_cd16_ifn_pic.obs['group'] = np.select(conditions, choices, default='Resting')


# Load top De genes in each state as a list
my_file = open(r"C:\Users\TzachiHNB6\Documents\annotated_count_matrices\bat\blood\bat1k_DE\cd16_mafb_de\bat_no_orthologs_de\bat_top_100_gene_groups_python.txt", "r")
data = my_file.readlines()
data = [gene.replace('\n', "") for gene in data]

# Set colors for cell gropus
groups = bat1k_all_data_cd16_ifn_pic.obs[["group"]]
lut = dict(zip(groups["group"].unique(), ["#c2c2d6","#bf80ff","#4d94ff"]))
col_colors = groups["group"].map(lut)

# Load gene-state table
gene_groups = pd.read_csv(r"C:\Users\TzachiHNB6\Documents\annotated_count_matrices\bat\blood\bat1k_DE\cd16_mafb_de\bat_no_orthologs_de\gene_groups_map.csv")
gene_groups.set_index("gene", inplace=True)


# Set colors for gene gropus
lut2 = dict(zip(gene_groups["group"].unique(), ["#c2c2d6","#4d94ff","#bf80ff"]))
row_colors = gene_groups["group"].map(lut2)

# Plot genes and cells using clustering 
# This is a non-deterministic mathod. clustering results may differ based on order of gene passed to method.
cg = sns.clustermap(bat1k_all_data_cd16_ifn_pic.to_df()[data].T,cmap = "mako", row_cluster=True,col_cluster=True, method="ward",vmin = -3 , vmax=3, col_colors=col_colors, row_colors=row_colors, xticklabels=False)
cg.ax_row_dendrogram.set_visible(False)
cg.ax_col_dendrogram.set_visible(False)
cg.ax_heatmap.set_ylabel("")

# Plot and save to image
plt.show()
plt.savefig(r"C:\Users\TzachiHNB6\Documents\annotated_count_matrices\bat\blood\bat1k_DE\cd16_mafb_de\bat_no_orthologs_de\heatmap.pdf")


