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
import sys
from cell2location.utils import select_slide

results_folder = './results2/'
run_name = f'{results_folder}/cell2location_map_adataa_lower_left/'

adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select up to 6 clusters
clust_labels = ['Stem_cell', 'TA', 'Enterocyte', 'SMC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adataa_lower_left.png")
plt.savefig("bat_adataa_lower_left.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatac_right'


adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatac_right.png")
plt.savefig("bat_adatac_right.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adataa_lower_right'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adataa_lower_right.png")
plt.savefig("bat_adataa_lower_right.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adataa_upper_left'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adataa_upper_left.png")
plt.savefig("bat_adataa_upper_left.pdf", format="pdf", bbox_inches="tight")


run_name = f'{results_folder}/cell2location_map_adataa_upper_right'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adataa_upper_right.png")
plt.savefig("bat_adataa_upper_right.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatab_middle'

adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatab_middle.png")
plt.savefig("bat_adatab_middle.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatab_upper'

adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatab_upper.png")
plt.savefig("bat_adatab_upper.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatac_left'


adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatac_left.png")
plt.savefig("bat_adatac_left.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatad_left'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatad_left.png")
plt.savefig("bat_adatad_left.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatad_lower'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatad_lower.png")
plt.savefig("bat_adatad_lower.pdf", format="pdf", bbox_inches="tight")

run_name = f'{results_folder}/cell2location_map_adatad_upper'
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
# select up to 6 clusters


with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=adata_vis,
        color=clust_col, labels=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right',
        reorder_cmap=(3,1,0,2)
        # Use custom color map
    )
plt.savefig("bat_adatad_upper.png")
plt.savefig("bat_adatad_upper.pdf", format="pdf", bbox_inches="tight")
 
