import os

input_folder = "/tzachi_storage/yanl/bat_genes/Masked_MSA/Phy/"
output_folder = "/tzachi_storage/yanl/bat_genes/Masked_MSA/log/"

protein_names = ["C_2", "C_3", "C_5", "C1S", "C1R", "C8g", "C8b", "C8a", "CFB", "CFD", "CFI", "C_9", "C_6", "C_7"]

for protein_name in protein_names:
    seqfile = f"{input_folder}{protein_name}.phy"
    treefile = f"{input_folder}{protein_name}.phy_phyml_tree_edited.txt"
    outfile = f"{output_folder}{protein_name}_MSA_positive.codeml"
    
    content = f"""seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 9
verbose = 1
runmode = 0
Mgene = 0
seqtype = 1
clock = 0
model = 0
NSsites = 8
ncatG = 6
icode = 0
aaDist = 0
CodonFreq = 2
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0
Malpha = 0
fix_alpha = 1
alpha = 0
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
fix_blength = 0
method = 1"""

    output_file_name = f"{output_folder}{protein_name}_MSA_positive.ctl"
    with open(output_file_name, "w") as f:
        f.write(content)

input_folder = "/tzachi_storage/yanl/bat_genes/Masked_MSA/Phy/"
output_folder = "/tzachi_storage/yanl/bat_genes/Masked_MSA/log/"

for protein_name in protein_names:
    seqfile = f"{input_folder}{protein_name}.phy"
    treefile = f"{input_folder}{protein_name}.phy_phyml_tree_edited.txt"
    outfile = f"{output_folder}{protein_name}_MSA_no_positive.codeml"
    
    content = f"""seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 9
verbose = 1
runmode = 0
Mgene = 0
seqtype = 1
clock = 0
model = 0
NSsites = 8
ncatG = 6
icode = 0
aaDist = 0
CodonFreq = 2
fix_kappa = 0
kappa = 2
fix_omega = 1
omega = 1
Malpha = 0
fix_alpha = 1
alpha = 0
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
fix_blength = 0
method = 1"""

    output_file_name = f"{output_folder}{protein_name}_MSA_no_positive.ctl"
    with open(output_file_name, "w") as f:
        f.write(content)
