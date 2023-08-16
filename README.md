Spatial and single-cell transcriptomics illuminate bat immunity and barrier tissue evolution
============================================================================================

This directory contains data files and scripts related to our work: 
"Spatial and single-cell transcriptomics illuminate bat immunity and barrier tissue evolution".

The scripts and documentations in this directory allow reproducing the processing of raw sc RNAseq ,the statistical analysis and the main plots. 

1) sc-rna-seq folder contains scripts for processing RNAseq using 10X CellRanger pipeline, as well as spatial data processing using 10X SpaceRanger and 
trajectory analysis using velocyto.

2) EggNOG folder contains scripts for runnig EggNOG annotations mapper and processing its results for producing orthology relations between species. 

3) cell2location folder contains scripts for using cell2location for cell type mapping in spatial transcriptomics.

4) cross-species folder contains scripts for integrating cross-species datasets and performing correlation analysis of them.

5) differential_expression folder contains scripts for conducting DE analyses for both signle-cell and bulk data.

6) coding_sequence_evolution folder contains a pipeline for running positive selection analysis on multi-sequense alighnment data.

7) figures folder contains scripts for reproducing the main plots.