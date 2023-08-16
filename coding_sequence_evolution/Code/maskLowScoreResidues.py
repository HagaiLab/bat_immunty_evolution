import sys
import os
import pbs_runners_clust
import file_utilities
from os import listdir
from os.path import isfile, join
import subprocess
# eggnog command 

msa_directory = "/tzachi_storage/yanl/bat_genes/MSA/"
output_directory = "/tzachi_storage/yanl/bat_genes/Masked_MSA"
os.makedirs(output_directory, exist_ok=True)
folder_list = sorted(os.listdir(msa_directory))

for folder_name in folder_list[:40]:  
  if "MSA" in folder_name:
    continue  

  folder_path = os.path.join(msa_directory, folder_name)
  msa_file = os.path.join(msa_directory, folder_name, "MSA.MAFFT.Guidance2_res_pair_res.scr")
  output_file = os.path.join(output_directory, folder_name + ".txt")

  cmd = f"""module load gcc/gcc-12.1.0
  module load perl/perl-5.30.1-new
  module load guidance/guidance-2.02-new
  module load mafft/7.450
  #cd /powerapps/share/guidance.v2.02/www/Guidance/
  cd /tzachi_storage/yanl/bat_genes/
  # Update the file paths based on the folder name
  perl maskLowScoreResidues.pl {os.path.join(msa_directory, folder_name, 'MSA.MAFFT.aln.With_Names')} {msa_file} {output_file} 0.9 nuc"""
  jobID = pbs_runners_clust.script_runner(cmds=cmd, alias="cds", load_python=True,
                                            gmem=2, ncpus=1, queue="adistzachi")
  print(jobID)
