
library(Seurat)
library(SeuratDisk)

#load human data - convert from Python H5AD object to Seurat Object
Convert("C:/Users/TzachiHNB6/Documents/annotated_count_matrices/human/gut/before_regress/new_human_after_ort.h5ad", dest = "C:/Users/TzachiHNB6/Documents/annotated_count_matrices/human/gut/before_regress/new_human_after_ort.h5seurat", overwrite = TRUE)
human_gut_data = LoadH5Seurat("C:/Users/TzachiHNB6/Documents/annotated_count_matrices/human/gut/before_regress/new_human_after_ort.h5seurat")
human_gut_data_copy <- human_gut_data

#load bat data - convert from Python H5AD object to Seurat Object
Convert("C:/Users/TzachiHNB6/Documents/annotated_count_matrices/bat/gut/before_regress/bat_50_after_ort.h5ad", dest = "C:/Users/TzachiHNB6/Documents/annotated_count_matrices/bat/gut/before_regress/bat_50_after_ort.h5seurat" ,overwrite = TRUE)
bat_gut_data <- LoadH5Seurat("C:/Users/TzachiHNB6/Documents/annotated_count_matrices/bat/gut/before_regress/bat_50_after_ort.h5seurat")
bat_gut_data_copy <- bat_gut_data

#Log-normalize and find High Variable Genes 
bat_gut_data_copy <- NormalizeData(bat_gut_data_copy, normalization.method = "LogNormalize", scale.factor = 10000)
bat_gut_data_copy <- FindVariableFeatures(bat_gut_data_copy, selection.method = "vst", nfeatures = 2000)

human_gut_data_copy <- NormalizeData(human_gut_data_copy, normalization.method = "LogNormalize", scale.factor = 10000)
human_gut_data_copy <- FindVariableFeatures(human_gut_data_copy, selection.method = "vst", nfeatures = 2000)

#set type to be species
bat_gut_data_copy$type <-"bat"
human_gut_data_copy$type <-"human"

# Integration
objects.list <- c(bat_gut_data_copy, human_gut_data_copy)
data.anchors <- FindIntegrationAnchors(object.list = objects.list)
gc()
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:30)
data.integrated <- ScaleData(data.integrated)
data.integrated <- RunPCA(data.integrated, npcs = 30)
data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 0.5)
gc()

#Save integrated object 
saveRDS(data.integrated, file = "annotated_count_matrices/seurat/data.rds")

#plot clusters by type and cell_type
DimPlot(data.integrated, reduction = "umap", group.by = c('type'))
DimPlot(data.integrated, reduction = "umap", group.by = c('inter_species_ann_B'), label =T)
