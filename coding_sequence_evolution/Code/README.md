Order of pipeline
============================================================================================
some of the steps are done using bash scripts and some using Python scripts

1) Create multi-sequence file for each gene

2) Run GUIDANCE for multi-sequence-alignment (guidance.txt)

3) Optional: run masking of residues (mask_low_score_residues.txt)

4) Convert MSA file  - masked or unmasked - to phylip (convert_fasta_to_phylip.py)
	
5) Optional: create phylogenetic tree(s) with phyml (tree_with_phyml.txt) / use your own tree(s)

6) Adjust tree for PAML use (tree_remove_numbers.py)

7) Create .ctl configuration file(s) for PAML (create_codeml_ctl_files.py)
	
8) run PAML's codeml (run_codeml.py)
	
