import sys
import os
import pbs_runners_clust
import file_utilities
from os import listdir
from os.path import isfile, join
import subprocess

_path = "/tzachi_storage/yanl/bat_genes/MAFFT_phyml" 

for file_name in [

'C1R_MSA_no_positive.ctl',
'C1R_MSA_positive.ctl',
'C1S_MSA_positive.ctl',
'C1S_MSA_no_positive.ctl',
'c8g_MSA_no_positive.ctl',
'c8g_MSA_positive.ctl',
'c8b_MSA_positive.ctl',
'c8b_MSA_no_positive.ctl',
'c8a_MSA_positive.ctl',
'c8a_MSA_no_positive.ctl',
'CFD_MSA_positive.ctl',
'CFD_MSA_no_positive.ctl',
'CFI_MSA_positive.ctl',
'CFI_MSA_no_positive.ctl',
'CFB_MSA_positive.ctl',
'CFB_MSA_no_positive.ctl',
'C_5_MSA_positive.ctl',
'c_7_MSA_positive.ctl',
'C_5_MSA_no_positive.ctl',
'c_7_MSA_no_positive.ctl',
'C_2_MSA_positive.ctl',
'C_2_MSA_no_positive.ctl',
'C_6_positive.ctl',
'C_6_no_positive.ctl',
'C_3_MSA_no_positive.ctl',
'C_3_MSA_positive.ctl',
'C_9_no_positive.ctl',
'C_9_positive.ctl'

]:

  cmd = """module load miniconda/miniconda3-4.7.12-environmentally 
  module load paml/paml4.9h 
  codeml {}""".format(join(_path, file_name))
  jobID = pbs_runners_clust.script_runner(cmds=cmd, alias="s_{}".format(file_name)[0:22], load_python=True, gmem=2, ncpus=1, queue="adistzachi") 
  print(jobID) 
