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
from cell2location import run_colocation

data_type = 'float32'
#module load cell2location

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')


print('loading done')

#rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
#print('font choosen')
results_folder = './results2/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map_adataa_upper_left'

adata_ref = sc.read('/tzachi_storage/yanl/Spatial_Cell2Loc/gut_raw_bat1k.h5ad')
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type'
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        #categorical_covariate_keys=['Method']
                       )
# create the regression model

mod = RegressionModel(adata_ref)
mod.train(max_epochs=250, use_gpu=True)
mod.plot_history(20)
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
 
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
print(inf_aver.iloc[0:5, 0:5])

#Cell2location: spatial mapping
#find shared genes and subset both anndata and reference signatures
adata_vis = sc.read('/tzachi_storage/yanl/Spatial_Cell2Loc/adataa_upper_left.h5ad')
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)#, batch_key="in_tissue")
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=18,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
#mod.view_anndata_setup()
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)


#mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
#fig = mod.plot_spatial_QC_across_batches()
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(12, 15), # IMPORTANT: use a wider range of the number of factors (5-30)
      #'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': f'{run_name}/CoLocatedComb/'}
)

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata_vis, cmap='magma',
                  # show first 8 cell types
                  color=['B', 'B_Plasma', 'Basophil/Eosinophil', 'Cycling_B', 'Cycling_B_Plasma', 'Cycling_T', 'ECs_cap/vas', 'EECs', 'Enterocyte', 'Enterocyte_HLA-A', 'Fibroblast_ADAMDEC1', 'Myofibroblast', 'Fibroblast_PI16', 'Glia', 'Goblet', 'ILC2', 'ILC3', 'LEC', 'Macrophage', 'Mast', 'Mesothelium', 'Monocyte', 'NK', 'NKT', 'Neutrophil', 'Paneth', 'Pericyte', 'SMC', 'Stem_cell', 'TA', 'T_CD4', 'T_CD8', 'T_CD8_GZMK', 'Tuft', 'cDC1', 'cDC2', 'pDC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
