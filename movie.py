"""
@author: Joe Laforet
Created on Sat June 21 14:41:50 2022


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""


import pandas as pd
import numpy as np
import os
import glob
import argparse
#####################################################

def writeJob(jobName, command):
    
    
    
    job_file = os.path.join(os.getcwd(), '%s.sh' % str(jobName))
    
    with open(job_file, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write('#SBATCH --mail-type=begin\n')
        fh.write('#SBATCH --mail-type=end\n')
        fh.write('#SBATCH --mail-user=jrl78@duke.edu\n')
        fh.write('#SBATCH -e slurm.err\n')
        fh.write("#SBATCH --job-name=%s\n" % jobName)
        fh.write("#SBATCH --mem=20G\n")
        fh.write("#SBATCH -p scavenger-gpu --gres=gpu:1\n")
        fh.write("#SBATCH --exclusive\n")
        
        fh.write('%s\n' % command)
	fh.write('find . -name "*.png" -type f -delete')


###################################################
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--indir", type=str, help="Input name of folder containing images.", required=True)
parser.add_argument("-d", "--del", action="store_true", default = False, help="If activated, will delete images after rendering video.")
parser.add_argument('-n', '--name', type=str, default = 'Test', help = 'Name for movie.', required = True)
parser.add_argument('-f', '--framerate', type=int, help='Framerate for output movie.', required = True)
args = parser.parse_args()
########################################################

framerate = args.framerate
name = args.name
os.chdir(name)
command = "ffmpeg -r {} -i {}_%04d.png -c:v libx264 -pix_fmt yuv420p -y {}.mp4".format(framerate, name, name)

writeJob(name, command)

os.system("chmod +x %s.sh" % str(name))
os.system("sbatch %s.sh" % str(name))

print("Job's Done!")

