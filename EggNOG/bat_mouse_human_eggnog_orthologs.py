import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
from Bio import SeqIO as si

"""
 
This script creates a list of 1:1:1 orthologs gene between human, mouse and bat (Rousettus aegyptiacus).
The script takes CDS fasta file and EggNOG's run results for each species 
and returns a dataframe containin all genes that have an identical EggNOG gene name. 
 
"""


def extract_cds_data(cds_path, gene_id_regex = None, gene_name_regex = None):   
    """
    Iterates CDS records and extracts record id (that will be used for intersection with Eggnog annotations file),
    gene id and gene name (if exists) 
    @param cds_path: cds file path
    @param gene_id_regex: pattern for extracting gene id from CDS record (optional)
    @param gene_name_regex: pattern for extracting gene name from CDS record (optional)
    @return: dataframe containing genes information
    """
    
    gene_ids = []
    gene_names = []
    fasta_record_ids = []
    with open(cds_path, "r") as cds:
        sequences = si.parse(cds, "fasta")
        for fasta in sequences:
            fasta_record_ids.append(fasta.id)
            if gene_id_regex != None:
                gene_ids.append(re.search(gene_id_regex, fasta.description).group(1))
            if gene_name_regex != None:
                gene_names.append(re.search(gene_name_regex, fasta.description).group(1))
    if gene_id_regex != None and gene_name_regex != None:
        return pd.DataFrame({'record_id' : fasta_record_ids, 'gene_id': gene_ids, 'gene_name': gene_names})
    if gene_id_regex != None:
        return pd.DataFrame({'record_id' : fasta_record_ids, 'gene_id': gene_ids})
    return pd.DataFrame({'record_id' : fasta_record_ids, 'gene_name': gene_names})
    


def combine_with_eggnog(data,annotations, gene_id_col_name = "record_id"):
    """
    Adds EggNOG gene name (if exist) for each gene
    @param data: dataframe containing genes information
    @param annotations: path to EggNOG annotations result file
    @param gene_id_col_name: column name in data represting record/gene id for merging with data from annoation file
    @return: updated dataframe containing EggNOG gene names
    """

    for index,row in data.iterrows():
        values = annotations.loc[annotations["query"] == row[gene_id_col_name], "Preferred_name"].values
        if len(values) > 0:
            data.loc[index,"eggnog_name"] =  values[0]
    return data


def clean_data(data):
    """
    Removes genes that have no EggNOG name and convert lower-case gene ids and names to upper-case
    @param data: dataframe containing genes information
    @return: updated dataframe
    """
    data = data[(data["eggnog_name"].notnull())& (data["eggnog_name"] != "-")]
    data["eggnog_name"] = [x.upper() for x in data["eggnog_name"]]
    data["gene_name"] = [x.upper() for x in data["gene_name"]]
    data.drop(columns = ["record_id"], inplace=True) # no need for that column anymore
    return data


def remove_multiple_gene_ids(data, gene_id_col = "gene_id"):
    """
    remove any gene that has multiple gene ids 
    @param data: dataframe containing genes information
    @return: updated dataframe containing only genes that have a single gene_id
    """
    data_grouped = data.groupby("eggnog_name").agg(set).reset_index() # group by eggnog name
    data_grouped = data_grouped[data_grouped[gene_id_col].map(len) == 1] # remove genes with more than 1 eggnog name
    data_grouped[gene_id_col] = [list(x)[0] for x in data_grouped[gene_id_col]]
    return data_grouped



def run():

    # Load annotation files from Eggnog run results.
    mouse_annotations = pd.read_csv(r"C:\Users\TzachiHNB6\Documents\eggnog\mouse_query_92_1.fa.emapper.annotations", delimiter="\t")
    human_annotations  = pd.read_csv(r"C:\Users\TzachiHNB6\Documents\eggnog\human_query_92_1.fa.emapper.annotations", delimiter="\t")
    bat1k_annotations = pd.read_csv(r"C:\Users\TzachiHNB6\Documents\eggnog\bat1k_query_1.fa.emapper.annotations", delimiter="\t")
    
    # Load CDS files and extract gene data
    cds_path = r"C:\Users\TzachiHNB6\Documents\eggnog\genoms\mouse_92_cds\Mus_musculus.GRCm38.cds.all.fa"
    mouse_data = extract_cds_data(cds_path, "gene:([^ ]+)", "gene_symbol:([^ ]+)")
    cds_path = r"C:\Users\TzachiHNB6\Documents\eggnog\genoms\human_92_cds\Homo_sapiens.GRCh38.cds.all.fa"
    human_data = extract_cds_data(cds_path, "gene:([^ ]+)", "gene_symbol:([^ ]+)")
    cds_path = r"C:\Users\TzachiHNB6\Documents\eggnog\genoms\bat_cds\cds_from_genomic.fna"
    bat_data = extract_cds_data(cds_path, gene_name_regex = "gene=([^\]]+)]")
    
    # Add eggnog annotations to gene data
    bat_data = combine_with_eggnog(bat_data, bat1k_annotations)
    human_data = combine_with_eggnog(human_data, human_annotations)
    mouse_data = combine_with_eggnog(mouse_data, mouse_annotations)

    bat_data_copy = bat_data.copy()
    human_data_copy = human_data.copy()
    mouse_data_copy = mouse_data.copy()
    
    # Clean data
    bat_data_copy = clean_data(bat_data_copy)
    mouse_data_copy = clean_data(mouse_data_copy)
    human_data_copy = clean_data(human_data_copy)

    # Remove version from gene id in human and mouse
    mouse_data_copy["gene_id"] = [x.split(".")[0] for x in mouse_data_copy["gene_id"]]
    human_data_copy["gene_id"] = [x.split(".")[0] for x in human_data_copy["gene_id"]]

    # Remove_multiple_gene_ids(bat_data_copy,"gene_name")
    bat_data_copy = remove_multiple_gene_ids(bat_data_copy,"gene_name")
    human_data_copy = remove_multiple_gene_ids(human_data_copy,"gene_id")
    mouse_data_copy  = remove_multiple_gene_ids(mouse_data_copy,"gene_id")

    # Change the gene name from a single string set to a single string
    human_data_copy["gene_name"] = [list(x)[0] for x in human_data_copy["gene_name"]]
    mouse_data_copy["gene_name"] = [list(x)[0] for x in mouse_data_copy["gene_name"]]


    # Change gene_name and gene_id column names to be species specific
    bat_data_copy.rename(columns={"gene_name": "bat_gene_name"}, inplace=True)
    human_data_copy.rename(columns={"gene_name": "human_gene_name", "gene_id":"human_gene_id"}, inplace=True)
    mouse_data_copy.rename(columns={"gene_name": "mouse_gene_name", "gene_id":"mouse_gene_id"}, inplace=True)

    merged = pd.merge(pd.merge(bat_data_copy,mouse_data_copy,on='eggnog_name', how="inner"),human_data_copy,on='eggnog_name', how="inner")     # Merge 3 genes datasets mkeeping only genes that contain an eggnog_name in all species, to get only genes that are 1:1:1 orthologs according to EggNOG DB

    merged.to_csv(r"bat_mouse_human_integration\bat_mouse_human_orthologs.csv") # Save 1:1:1 orthologs to file

if __name__ == "__main__":
    run()