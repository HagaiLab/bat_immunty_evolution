library(Seurat)
library(SeuratDisk)

#Here we perform DE-analysis between bat and human gut data in each of their shared cell types.

#Load integrated bat-human gut data
data.integrated = readRDS("bat_human_integrated_gut.rds")

#Set the identity to be the species 
Idents(data.integrated) = data.integrated$type

#Path of the folder for saving results.
base_path = "bat/gut/integration/DE/"

for(cell_type in unique(data.integrated$inter_species_ann)){
  # get the specific cell type data from the integrated data
  cell_type_subset = subset(data.integrated, subset = (inter_species_ann == cell_type))
  de_results = DE(cell_type_subset,"bat", "human")
  cell_type_corrected = gsub("/","_",cell_type) #replace "/" in cell type
  file_path = paste0(base_path ,sprintf("%s.csv",cell_type_corrected))
  write.csv(de.markers, file =file_path)
 
}

DE = function(data, ident.1, ident.2){
   de.markers <- FindMarkers(data ,ident = ident.1 , ident.2 = ident.2, assay ="RNA")
  #de.markers <- FindMarkers(data, min.pct = 0 ,logfc.threshold = 0 ,ident = ident.1 , ident.2 = ident.2, assay ="RNA") # No defalut restrictions
   return(de.markers)
}