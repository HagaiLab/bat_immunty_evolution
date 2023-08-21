#!/bin/bash
# here we use PAML for running phylogenetic analyses of sequences using maximum likelihood.
# PAML is used after creating the .ctl configuration files.
# installation, Documentation and Usage of PAML can be found in: http://abacus.gene.ucl.ac.uk/software/paml.html

CTL_Path="/tzachi_storage/yanl/bat_genes/MAFFT_phyml"

for entry in "$CTL_Path"/*
do
    codeml $entry
done


