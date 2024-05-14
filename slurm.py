import pandas as pd
import numpy as np
import os
import glob
import argparse


def writeJob(jobName, command):
    
    
    
    job_file = os.path.join(os.getcwd(), '%s.sh' % str(jobName))
    
    with open(job_file, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write('#SBATCH --mail-type=begin\n')
        fh.write('#SBATCH --mail-type=end\n')
        fh.write('#SBATCH --mail-user=jrl78@duke.edu\n')
        fh.write('#SBATCH -e slurm.err\n')
        fh.write("#SBATCH --job-name=%s.job\n" % jobName)
        fh.write("#SBATCH --mem=2G\n")
        fh.write("#SBATCH -p scavenger-gpu --gres=gpu:1\n")
        fh.write("#SBATCH --exclusive\n")
        
        fh.write('%s\n' % command)
       

    
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=str, help="Input command to be executed.", required=True)
parser.add_argument("-r", "--run", action="store_true", default = False, help="If activated, will send job script to cluster via sbatch.")
parser.add_argument('-n', '--name', type=str, default = 'Test', help = 'Name of job file.', required = True)
args = parser.parse_args()    

command = args.infile
name = args.name

writeJob(name, command)

print("Job written!")

if args.run == True:
    os.system("sbatch %s.sh" % str(name))
    print("Job sent!")