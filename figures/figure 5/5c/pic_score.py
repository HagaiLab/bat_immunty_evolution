import scanpy as sc
import anndata as ad
import pandas as pd

"""

 Primary innate response score calculation in R.aegyptiacus PBMC's
 Scoring is done on Poly(I:C) and LPS infected cells within chosen cell types
 Same scoring method was used in other scores analyses

"""

def load_h5ad(adata_path, metadata_path=None):
    """
    Load anndata object and metadata
    @param adata: anndata file path
    @param metadata_path: obs metadata (optional)
    @return adata: anndata object
    """
    adata = sc.read(adata_path)
    if metadata_path != None:
        bat_metadata = pd.read_csv(metadata_path) # Load bat PBMC's metadata
        adata.obs  = bat_metadata
    return adata

def clean_data(adata):
    """
    Filter out non releveant cells
    @param: anndata object
    @return: filtered anndata object
    """
    adata = adata[adata.obs.inter_species_ann.isin([ 'Neutrophils',
    'Monocytes_CD14',
    'Monocytes_CD16',
    'Monocytes_IFNB1',
    'DC_activated',
    'pDCs',
    'T_CD4',
    'T_CD8/NKT',
    'NK',
    'B',
    'B_CD83',
    ])]

    # Keep only LPS/Poly(I:C) infected cells
    adata = adata[adata.obs["batch"].isin(["LPS","Poly(I:C)"])]
    return adata

    
    
def parse_de(file_path):
    """
    Load DE results between Poly(I:C) and CTRL cells done in different species
    @param adata: DE file path
    @return de_pic_ctrl: dataframe of DE results
    """
    de_pic_ctrl = pd.read_excel(file_path)
    # Only keep human, mouse and R.aegyptiacus) DE results
    de_pic_ctrl = de_pic_ctrl[["Homo sapiens_ortho_gene","Mus musculus_ortho_gene","Rousettus aegyptiacus_ortho_gene","Homo sapiens FC","Homo sapiens Q-value","Mus musculus FC","Mus musculus Q-value","Rousettus aegyptiacus FC","Rousettus aegyptiacus Q-value"]]
    return de_pic_ctrl
  
  
def get_scoring_gene_list(de_data):
    """
    Keeps only genes that are upregulated (significantly) in Poly(I:C) compared with CTRL in releveant species
    @param de_data:  dataframe of DE results
    @return: a list of genes to score by
    """
    pic_genes = []
    # Pair up FC and Q-values of each relevant species
    fc_qval_pairs = {"Homo sapiens FC":"Homo sapiens Q-value", "Rousettus aegyptiacus FC":"Rousettus aegyptiacus Q-value", "Mus musculus FC": "Mus musculus Q-value" }
    for index, row in de_data.iterrows():
        passed = True
        for fc, qval in fc_qval_pairs.items():       
            if (row[fc] <0) | (row[qval]  > 0.05):
                passed = False
                break
        if passed:
            bat_gene = row.loc["Rousettus aegyptiacus_ortho_gene"].upper()
            pic_genes.append(bat_gene)
    return pic_genes
            
            
def score_genes(adata,genes,score_title):
    """
    Score calculating function (using sc.tl.score_genes)
    @param adata: anndata object
    @param genes: a list of genes to calculate score on
    @param score_title: string title for plotting
    @return adata: anndata object with scores in obs
    """
    adata_copy = adata.copy()
    sc.pp.normalize_per_cell(adata_copy, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.scale(adata_copy)
    sc.tl.score_genes(adata_copy,gene_list=genes, score_name= score_title+'_score')
    adata_cc_genes = adata_copy[:, genes].copy()
    adata.obs[score_title+' score'] = adata_copy.obs[score_title+'_score'].copy()
    return adata
    
    
def run():
    # Load bat PBMC's data
    adata = load_h5ad(r"C:\Users\TzachiHNB6\Documents\annotated_count_matrices\bat\blood\bat_mouse_integration\PBMCs_annotated_raw_bat1k\PBMCs_annotated_raw_bat1k.h5ad",r"C:\Users\TzachiHNB6\Downloads\pbmcs_bat1k_metadata.csv")
    adata_copy = adata.copy()
    
    #Clean data
    adata_copy = clean_data(adata_copy)
    
    #Load and parse de data
    de_data = parse_de(r"C:\Users\TzachiHNB6\Downloads\pIC_up_six_mammals_fiborblasts.xlsx")
   
    # Extract genes for scoring
    pic_genes = get_scoring_gene_list(de_data)
    
    # Find genes names in list mismatching names in R.aegyptiacus genes and rename them
    adata_copy.var_names = [i.upper() for i in adata_copy.var_names]
    print([i for i in pic_genes if i not in adata_copy.var_names]) 
    replace = {"HERC6":"HERC6_LOC107518925", "LOC107519342":"ZC3HAV1", "FAM46A":"TENT5A", "GSAP":"GSAP_LOC107511514","MB21D1" :"CGAS", 
               'LOC107500491':"SPACA6", "LOC107515132":"SHFL"}  # Correct mismatching gene names in gene list
    pic_genes = [replace[x]  if x in replace else x for x in pic_genes]
    pic_genes.remove("IRF1")
    # Make sure no gene name is mismatched
    print([i for i in pic_genes if i not in adata_copy.var_names])

    # Run scoring
    score_genes(adata_copy,pic_genes,'PIC')

    # Plot scores and save to pdf (craetes a "figures" folder in working directory and saves pdf in it) 
    sc.pl.matrixplot(adata_copy ,'PIC score',groupby ='inter_species_ann',categories_order = ['Neutrophils','Monocytes_CD14','Monocytes_CD16','Monocytes_IFNB1','DC_activated','pDCs','T_CD4','T_CD8/NKT','NK','B','B_CD83'],dendrogram=False,standard_scale='var',cmap='coolwarm',colorbar_title='Scaled \nscore',save='bat_pic_score_plot.pdf')


if __name__ == "__main__":
    run()