#! /usr/local/python_anaconda/bin/python3.4

import pbs_jobs
import os
from os import path
from file_utilities import set_filenames_for_pbs_runs, check_filename, check_dirname
import pandas as pd


def script_runner(cmds, alias = "Bulk_quant_21_03_21", load_python=True, gmem=2, ncpus=1, queue="adistzachi"):
    """
    run script on cluster
    :param cmds: script running line
    :param alias: job name (default: script)
    :return: job id
    """
    cmdfile = pbs_jobs.get_cmdfile_dir("scriptBulk_quant_21_03_21", alias)
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds, ncpus=ncpus, load_python=load_python)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id
