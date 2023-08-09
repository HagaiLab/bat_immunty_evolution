import os
import subprocess
from time import sleep
import getpass
import datetime

def create_pbs_cmd(cmdfile, alias, queue="dudulight", gmem=2, ncpus=1, ngpus=1, cmds="", dir = "",
                   load_python=True, jnum=False, run_after_job=None):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q {}\n".format(queue))
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N {}\n".format(alias))

        if alias in cmdfile and datetime.datetime.today().strftime('%Y-%m') in cmdfile:
            o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
            o.write("#PBS -e %s\n" % "/".join(cmdfile.split("/")[:-1]))

        # running on GPUs 
        if queue == 'gpu':
            o.write("#PBS -l select=ngpus={}\n".format(ngpus))
        else:    
            o.write("#PBS -l select=ncpus={}:mem={}gb\n".format(ncpus,gmem))

        if jnum:
           if jnum != 1:
               o.write("#PBS -J 1-{}\n\n".format(jnum))
        if run_after_job != None:
            o.write("#PBS -W depend=afterok: {}\n\n".format(run_after_job))
    
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("date\n")
        o.write("hostname\n")
        #if load_python:
         #   o.write("module load python/anaconda_python-3.6.1\n")
        #o.write("module load salmon-0.13.1\n")
        
        # for eggnog
        o.write("module load miniconda/miniconda3-environmentally\n")
        o.write("conda activate /powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/eggnog-2.1.2\n")
        
        
        o.write("\n")
        #o.write("echo {}".format(cmds))
        o.write("\n")
        o.write(cmds)
        o.write("\n")
        o.write("date\n")
    o.close()


def submit(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]


def check_pbs(job_id):
    """
    :param job_id: The PBS job id
    :return: "Done!", when the job is done
    """
    status = "Running..."
    try:
        process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
        while process != "":
            process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
            sleep(0.05)
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status


def get_cmdfile_dir(cmdfile, alias):
    username = getpass.getuser()
    lab_users_dic = {"lilachs": "/tzachi_storage/lilachs/eggnog/logs_eggnog_11_04_21"}
    if username in lab_users_dic.keys():
        tmp_dir = lab_users_dic[username]
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        date = datetime.datetime.today().strftime('%Y-%m')
        tmp_dir = tmp_dir + "/%s" % date
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        tmp_dir = tmp_dir + "/%s" % alias
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        cmdfile = tmp_dir + "/" + cmdfile
    return cmdfile
