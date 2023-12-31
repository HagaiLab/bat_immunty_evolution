import random
import anndata
import scanpy as sc
import pandas as pd

def down_sampling2(bdatas): # list of adatas - order needs to be - [adata_bat, adata_mouse,adata_human]
    adatas = bdatas.copy()
    barcodes = []
    numbers = {}
    
    # Remove cells with less than 10 counts from each dataset
#     for i in range(len(adatas)):
#         adatas[i] = adatas[i][~adatas[i].obs_names.isin(adatas[i][adatas[i].obs['inter_species_ann'].isin(set(adatas[i].obs['inter_species_ann'].value_counts()[adatas[i].obs['inter_species_ann'].value_counts() < 10].index))].obs_names)]
        
    # Concatenate the datasets
    adatall = adatas[0].concatenate(adatas[1:], batch_categories=['bat', 'mouse','human'],index_unique=None,)

    x = set(adatall[adatall.obs['batch'] == 'bat'].obs['inter_species_ann'])
    y = set(adatall[adatall.obs['batch'] == 'mouse'].obs['inter_species_ann'])
    z = set(adatall[adatall.obs['batch'] == 'human'].obs['inter_species_ann'])
    cell_types = (x & y) | (y & z) | (z & x)

    for celltype in cell_types:
        numbers[celltype] = 0
        for i in range(len(adatas)):
            if celltype in set(adatas[i].obs['inter_species_ann']):
                numbers[celltype] += adatas[i][adatas[i].obs['inter_species_ann'] == celltype].shape[0]

    random.seed(1)
    for celltype in cell_types:
        if numbers[celltype] > 0:
            sample_size = min([adata[adata.obs['inter_species_ann'] == celltype].shape[0] for adata in adatas if celltype in set(adata.obs['inter_species_ann'])])
            for adata in adatas:
                if celltype in set(adata.obs['inter_species_ann']):
                    for i in random.sample(list(adata[adata.obs['inter_species_ann'] == celltype].obs_names), k=sample_size):
                        barcodes.append(i)
        
    return adatas[0][adatas[0].obs.index.isin(barcodes)], adatas[1][adatas[1].obs.index.isin(barcodes)], adatas[2][adatas[2].obs.index.isin(barcodes)]



#blood
# b = sc.read(r"C:\Users\TzachiHNB7\Downloads\Bat\PBMCs_annotated_raw_bat1k.h5ad")#Bat
# b.obs.inter_species_ann= b.obs.inter_species_ann.astype('str')
# b.obs.cell_type= b.obs.cell_type.astype('str')
# b.obs.inter_species_ann.replace({'DCs_LAMP3_FSCN1':'DC_activated'},inplace=True)
# b.obs.cell_type.replace({'DCs_LAMP3_FSCN1':'DC_activated'},inplace=True)

# for i in b.obs_names:
#     if b.obs.at[i,'inter_species_ann'] == 'T_CD8/NKT':
#         if b.obs.at[i,'cell_type'] in ['T_CD8_CM','T_CD8_EM']: ##Replace cell_type names if needed in another tissue/species
#             b.obs.at[i,'inter_species_ann'] = 'T_CD8'
#         else: 
#             b.obs.at[i,'inter_species_ann'] = 'NKT'
# b.obs.inter_species_ann_B= b.obs.inter_species_ann+' (B)'

# h = sc.read(r"C:\Users\TzachiHNB7\Downloads\Human\human_pbmcs_raw_annotated.h5ad") # Human
# h.obs.inter_species_ann= h.obs.inter_species_ann.astype('str')
# h.obs.inter_species_ann_H= h.obs.inter_species_ann_H.astype('str')
# for i in h.obs_names:
#     if h.obs.at[i,'annotation_broad'] == 'T CD8+':
#         h.obs.at[i,'inter_species_ann'] = 'T_CD8'
#     if h.obs.at[i,'annotation_broad'] == 'T g/d':
#         h.obs.at[i,'inter_species_ann'] = 'T_g/d'
#     if h.obs.at[i,'annotation_broad'] == 'MAIT':
#         h.obs.at[i,'inter_species_ann'] = 'MAIT'
# h.obs.inter_species_ann_H= h.obs.inter_species_ann+' (H)'

# m = sc.read(r"C:\Users\TzachiHNB7\Downloads\Mouse\PBMCs_mice_raw.h5ad") # Mouse
# mouse_obs = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Mouse\mouse_pbmcs_metadata.csv", index_col=0)
# pdcs = mouse_obs[mouse_obs.cell_type=='pDCs'].index
# DC_activated = mouse_obs[mouse_obs.cell_type=='DC_activated'].index
# m.obs.inter_species_ann= m.obs.inter_species_ann.astype('str')
# m.obs.cell_type= m.obs.cell_type.astype('str')
# m.obs.inter_species_ann_M= m.obs.inter_species_ann_M.astype('str')

# pdcs = mouse_obs[mouse_obs.cell_type=='pDCs'].index
# DC_activated = mouse_obs[mouse_obs.cell_type=='DC_activated'].index
# m.obs.inter_species_ann= m.obs.inter_species_ann.astype('str')
# m.obs.cell_type= m.obs.cell_type.astype('str')
# m.obs.inter_species_ann_M= m.obs.inter_species_ann_M.astype('str')
# for i in m.obs_names:
#     if i in pdcs:
#         m.obs.at[i,'inter_species_ann'] = 'pDCs'
#         m.obs.at[i,'cell_type'] = 'pDCs'
#     elif i in DC_activated:
#         m.obs.at[i,'inter_species_ann'] = 'DC_activated'
#         m.obs.at[i,'cell_type'] = 'DC_activated'
#     elif i in m[m.obs.inter_species_ann=='T_CD8/NKT'].obs_names:
#         if m.obs.at[i,'cell_type'] == 'T_CD8':
#             print('T_CD8')
#             m.obs.at[i,'inter_species_ann'] = 'T_CD8'
#         else:
#             m.obs.at[i,'inter_species_ann'] = 'NKT'
#     else:
#         pass
        
# m.obs.inter_species_ann_M= m.obs.inter_species_ann+' (M)'

# Lung
# h = sc.read(r"C:\Users\TzachiHNB7\Downloads\Human\Elo_lung_both_sexes.h5ad")
# h = sc.read(r"C:\Users\TzachiHNB7\Downloads\Human\Human_Male_Lungs_Elo_last_version.h5ad")
# b = sc.read(r"C:\Users\TzachiHNB7\Downloads\Bat\Lung_annotated_raw_bat1k.h5ad")
# m = sc.read(r"C:\Users\TzachiHNB7\Downloads\Mouse\mouse_lungs_raw.h5ad")
# mouse_meta = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Mouse\mouse_lungs_metadata.csv", index_col = 0)
# human_meta = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Human\human_M_lung_metadata.csv",index_col = 0)
# bat_meta = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Bat\Bat1k_lungs_metadata.csv",index_col = 0)
# Gut
h = sc.read(r"C:\Users\TzachiHNB7\Downloads\Human\human_small_intestine_raw.h5ad")
b = sc.read(r"C:\Users\TzachiHNB7\Downloads\Bat\gut_raw_bat1k.h5ad")
m = sc.read(r"C:\Users\TzachiHNB7\Downloads\Mouse\mouse_small_intestine_raw.h5ad")
# Intestine

# Roy said NO to downsample in 5.2 !!! but do intersection between the inter_species_ann and no MT threshold amen.
dict_replace = {'COX1':'MT-COX1','COX2':'MT-COX2','COX3':'MT-COX3','ND1':'MT-ND1','ND2':'MT-ND2',
                 'ND3':'MT-ND3','ND4':'MT-ND4','ND5':'MT-ND5','ND6':'MT-ND6','ND4L':'MT-ND4L','ATP6':'MT-ATP6','ATP8':'MT-ATP8',
                 'CYTB':'MT-CYTB'} 
h.var.rename(dict_replace, inplace = True)
h.var['MT'] = h.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(h, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)
#h= h[h.obs.pct_counts_MT < 20, :] # redundant
#b = b[b.obs.pct_counts_MT < 20, :] # redundant
#m = m[m.obs.pct_counts_MT < 20, :] # redundant
h.obs.drop(columns = ['n_genes_by_counts', 'total_counts', 'total_counts_MT'], inplace=True)
del h.obsm
m.var = m.var.reset_index().set_index('Ensembl')
b.var = b.var.reset_index().set_index('original_name')

#lung
h = h[h.obs.index.isin(human_meta.index)]
h.obs['inter_species_ann_H'] =human_meta['inter_species_ann_H']
h.obs['inter_species_ann'] = human_meta['inter_species_ann']
h.obs['cell_type'] = human_meta['cell_type']
#Lung
b = b[b.obs.index.isin(bat_meta.index)]
b.obs['inter_species_ann_B'] = bat_meta['inter_species_ann_B']
b.obs['inter_species_ann'] = bat_meta['inter_species_ann']
b.obs['cell_type'] = bat_meta['cell_type']
#Lung
m = m[m.obs.index.isin(mouse_meta.index)]
m.obs['inter_species_ann_M'] = mouse_meta['inter_species_ann_M']
m.obs['inter_species_ann'] = mouse_meta['inter_species_ann']
m.obs['cell_type'] = mouse_meta['cell_type']

ort = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Orthologs\mouse_human_bat_orthologs.csv", index_col='mouse_gene_id')
#ort_bh = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Orthologs\bat1k_human_orthologs.csv", index_col='human gene id') # for correlations
#ort_bm = pd.read_csv(r"C:\Users\TzachiHNB7\Downloads\Orthologs\bat1k_mouse_orthologs.csv", index_col='mouse gene id') # for correlations

# IN FULL Only (cross - species)
temp = m.var[m.var.index.isin(ort.index)]
m = m[:, m.var.index.isin(temp.index)]
m.var['ind'] = ort['eggnog_name']
m.var = m.var.reset_index().set_index('ind')
ort = ort.reset_index().set_index('bat gene name')
temp = b.var[b.var.index.isin(ort.index)]
b = b[:, b.var.index.isin(temp.index)]
b.var['ind'] = ort['eggnog_name']
b.var = b.var.reset_index().set_index('ind')
ort = ort.reset_index().set_index('human gene id')

#h.var = h.var.reset_index().set_index('Ensembl')
#h.var = h.var.reset_index().set_index('gene_ids-1') #lung
h.var = h.var.reset_index().set_index('gene_ids') #gut

# Twins (correlations)
temp = m.var[m.var.index.isin(ort_bm.index)]
mb = m[:, m.var.index.isin(temp.index)]
mb.var['ind'] = ort_bm['eggnog_name']
ort_bm = ort_bm.reset_index().set_index('bat gene name')
temp = b.var[b.var.index.isin(ort_bm.index)]
bm = b[:, b.var.index.isin(temp.index)]
bm.var['ind'] = ort_bm['eggnog_name']

temp2 = h.var[h.var.index.isin(ort_bh.index)]
hb = h[:, h.var.index.isin(temp2.index)]
hb.var['ind'] = ort_bh['eggnog_name']
ort_bh = ort_bh.reset_index().set_index('bat gene name')
temp2 = b.var[b.var.index.isin(ort_bh.index)]
bh = b[:, b.var.index.isin(temp2.index)]
bh.var['ind'] = ort_bh['eggnog_name']

# IN FULL Only (cross - species)
temp = h.var[h.var.index.isin(ort.index)]
h = h[:, h.var.index.isin(temp.index)]
h.var['ind'] = ort['eggnog_name']
h.var = h.var.reset_index().set_index('ind')
h.var = h.var.drop(columns = 'index')

# Blood only (ALL TOGETHER)
temp = m.var[m.var.index.isin(h.var.index)]
m = m[:, m.var.index.isin(temp.index)]
temp2 = b.var[b.var.index.isin(h.var.index)]
b = b[:, b.var.index.isin(temp2.index)]

# Blood human mouse bat mouse (correlations)
hb.var = hb.var.reset_index().set_index('ind')
bh.var = bh.var.reset_index().set_index('ind')
mb.var = mb.var.reset_index().set_index('ind')
bm.var = bm.var.reset_index().set_index('ind')
temp = bh.var[bh.var.index.isin(hb.var.index)]
bat = bh[:, bh.var.index.isin(temp.index)]
temp = hb.var[hb.var.index.isin(bh.var.index)]
human = hb[:, hb.var.index.isin(temp.index)]
temp = bm.var[bm.var.index.isin(mb.var.index)]
bmm = bm[:, bm.var.index.isin(temp.index)]
temp = mb.var[mb.var.index.isin(bm.var.index)]
mbb = mb[:, mb.var.index.isin(temp.index)]
from scipy.sparse import csr_matrix
import anndata
df = h.to_df()
matrix = csr_matrix(df.values)
h.X = matrix
b, m, h = down_sampling2([b, m, h])
 
anndata.AnnData.write_h5ad(h, 'C:/Users/TzachiHNB7/Downloads/Human/human_bat1k_mouse_orthologs_gut_downsampled.h5ad')
anndata.AnnData.write_h5ad(b, "C:/Users/TzachiHNB7/Downloads/Bat/bat1k_mouse_human_orthologs_gut_downsampled.h5ad") 
anndata.AnnData.write_h5ad(m, 'C:/Users/TzachiHNB7/Downloads/Mouse/mouse_bat1k_human_orthologs_gut_downsampled.h5ad') 

### Correlations intersect + save###
common_cell_types = set(bmm.obs["inter_species_ann"]).intersection(set(mbb.obs["inter_species_ann"]))
# Filter the anndata objects to include only the common cell types
mbb = mbb[mbb.obs["inter_species_ann"].isin(common_cell_types)]
bmm = bmm[bmm.obs["inter_species_ann"].isin(common_cell_types)]
common_cell_types2 = set(human.obs["inter_species_ann"]).intersection(set(bat.obs["inter_species_ann"]))
# Filter the anndata objects to include only the common cell types
human = human[human.obs["inter_species_ann"].isin(common_cell_types2)]
bat = bat[bat.obs["inter_species_ann"].isin(common_cell_types2)]

# common_cell_types = set(b.obs["inter_species_ann"]).intersection(set(m.obs["inter_species_ann"]))
# # Filter the anndata objects to include only the common cell types
# mb = m[m.obs["inter_species_ann"].isin(common_cell_types)]
# bm = b[b.obs["inter_species_ann"].isin(common_cell_types)]
# common_cell_types2 = set(h.obs["inter_species_ann"]).intersection(set(b.obs["inter_species_ann"]))
# # Filter the anndata objects to include only the common cell types
# hb = h[h.obs["inter_species_ann"].isin(common_cell_types2)]
# bh = b[b.obs["inter_species_ann"].isin(common_cell_types2)]

anndata.AnnData.write_h5ad(bat, r"C:\Users\TzachiHNB7\Downloads\Bat\bat1k_human_orthologs_lung_not_ds.h5ad")
anndata.AnnData.write_h5ad(human, r"C:\Users\TzachiHNB7\Downloads\Human\human_bat1k_orthologs_lung_not_ds.h5ad")
anndata.AnnData.write_h5ad(bmm, "C:/Users/TzachiHNB7/Downloads/Bat/bat1k_mouse_orthologs_lung_not_ds.h5ad") 
anndata.AnnData.write_h5ad(mbb, 'C:/Users/TzachiHNB7/Downloads/Mouse/mouse_bat1k_orthologs_lung_not_ds.h5ad') 
