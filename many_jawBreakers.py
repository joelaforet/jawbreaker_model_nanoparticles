import argparse
import os
import subprocess

##################################################################################
parser = argparse.ArgumentParser(description='Generate jawBreakers for simulation')

parser.add_argument('-d', '--drug', nargs='+', required = True, help='list of drug names')
parser.add_argument('-e', '--excipient', nargs='+', required = True, help='list of excipient names')
parser.add_argument('--dnum', type=int, default=100, help='number of drug molecules')
parser.add_argument('-b', '--buffer', type=int, default=12, help='buffer distance (A)')

args = parser.parse_args()

#########################################################################

# Define the writeJob function
def writeJob(jobName, command):
    job_file = os.path.join(os.getcwd(), '%s.sh' % str(jobName))
    with open(job_file, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write('#SBATCH --mail-type=begin\n')
        fh.write('#SBATCH --mail-type=end\n')
        fh.write('#SBATCH --mail-user=jrl78@duke.edu\n')
        fh.write('#SBATCH -e slurm.err\n')
        fh.write("#SBATCH --job-name=%s.job\n" % jobName)
        fh.write("#SBATCH --mem=5G\n")
        fh.write("#SBATCH -p scavenger-gpu --gres=gpu:1\n")
        fh.write("#SBATCH --exclusive\n")
        fh.write('%s\n' % command)
    # Submit the job using sbatch
    subprocess.run(['sbatch', job_file])


for drug in args.drug:
    for excipient in args.excipient:
        jawbreaker_command = f"python make_jawBreaker.py -d {drug} -e {excipient} --DMSO DMSO --dNum {args.dnum} --solvent_type 4 -b {args.buffer}"
        writeJob(f"{drug}_{excipient}", jawbreaker_command)
