import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns
from cell2location.utils import select_slide
from cell2location.plt import plot_spatial
import matplotlib.colors as mcolors


results_folder = './results2/'
run_name = f'{results_folder}/cell2location_map_adatac_right'


adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']


# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatac_right1.png")
plt.savefig("bat_adatac_right1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatac_right2.png")
plt.savefig("bat_adatac_right2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatac_right3.png")
plt.savefig("bat_adatac_right3.pdf", format="pdf", bbox_inches="tight")


run_name = f'{results_folder}/cell2location_map_adataa_upper_left'


adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']


# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_upper_left1.png")
plt.savefig("bat_adataa_upper_left1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adataa_upper_left2.png")
plt.savefig("bat_adataa_upper_left2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_upper_left3.png")
plt.savefig("bat_adataa_upper_left3.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adataa_lower_left'


adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']


# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_lower_left1.png")
plt.savefig("bat_adataa_lower_left1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adataa_lower_left2.png")
plt.savefig("bat_adataa_lower_left2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_lower_left3.png")
plt.savefig("bat_adataa_lower_left3.pdf", format="pdf", bbox_inches="tight")

sys.exit()
run_name = f'{results_folder}/cell2location_map_adataa_lower_right'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_lower_right1.png")
plt.savefig("bat_adataa_lower_right1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adataa_lower_right2.png")
plt.savefig("bat_adataa_lower_right2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adataa_lower_right3.png")
plt.savefig("bat_adataa_lower_right3.pdf", format="pdf", bbox_inches="tight")
                              
run_name = f'{results_folder}/cell2location_map_adataa_upper_right'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adataa_upper_right1.png")
plt.savefig("bat_adataa_upper_right1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adataa_upper_right2.png")
plt.savefig("bat_adataa_upper_right2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adataa_upper_right3.png")
plt.savefig("bat_adataa_upper_right3.pdf", format="pdf", bbox_inches="tight")                              

run_name = f'{results_folder}/cell2location_map_adatab_middle'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatab_middle1.png")
plt.savefig("bat_adatab_middle1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatab_middle2.png")
plt.savefig("bat_adatab_middle2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatab_middle3.png")
plt.savefig("bat_adatab_middle3.pdf", format="pdf", bbox_inches="tight")  

run_name = f'{results_folder}/cell2location_map_adatab_upper'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatab_upper1.png")
plt.savefig("bat_adatab_upper1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatab_upper2.png")
plt.savefig("bat_adatab_upper2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatab_upper3.png")
plt.savefig("bat_adatab_upper3.pdf", format="pdf", bbox_inches="tight")  
            
run_name = f'{results_folder}/cell2location_map_adatac_left'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatac_left1.png")
plt.savefig("bat_adatac_left1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatac_left2.png")
plt.savefig("bat_adatac_left2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatac_left3.png")
plt.savefig("bat_adatac_left3.pdf", format="pdf", bbox_inches="tight") 

run_name = f'{results_folder}/cell2location_map_adatad_left'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatad_left1.png")
plt.savefig("bat_adatad_left1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatad_left2.png")
plt.savefig("bat_adatad_left2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatad_left3.png")
plt.savefig("bat_adatad_left3.pdf", format="pdf", bbox_inches="tight") 

run_name = f'{results_folder}/cell2location_map_adatad_lower'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatad_lower1.png")
plt.savefig("bat_adatad_lower1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatad_lower2.png")
plt.savefig("bat_adatad_lower2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatad_lower3.png")
plt.savefig("bat_adatad_lower3.pdf", format="pdf", bbox_inches="tight") 

run_name = f'{results_folder}/cell2location_map_adatad_upper'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
plt.savefig("bat_adatad_upper1.png")
plt.savefig("bat_adatad_upper1.pdf", format="pdf", bbox_inches="tight")
                
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )   
plt.savefig("bat_adatad_upper2.png")
plt.savefig("bat_adatad_upper2.pdf", format="pdf", bbox_inches="tight")
                             
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
                 
plt.savefig("bat_adatad_upper3.png")
plt.savefig("bat_adatad_upper3.pdf", format="pdf", bbox_inches="tight") 

#res_dict['n_fact12']['mod'].plot_cell_type_loadings()
