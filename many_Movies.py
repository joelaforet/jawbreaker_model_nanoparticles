"""
many_Movies.py
By Joe Laforet Jr.
jrl78@duke.edu

This script serves as an easy way to render multiple simulations at once. NOTE The .mol2 files of each molecule (drug/excipient) must also be in the current directory.

usage: many_Movies.py [-h]  
      optional arguments:
	  -h, --help            show this help message and exit
      -r REPLICATES,
                Input number of replicates to run for each simulation
      -t TIME,
                Input time to run simulation for. (Nanoseconds)
      -f FRAMES,
                Input amount of frames to save.
"""

import glob
import os
import mdtraj as md
import warnings
import subprocess
warnings.filterwarnings("ignore")

print('--------------------------------------------------------------------')
print("     Welcome to many_Movies.py!      ")
print("     By Joe Laforet Jr.      ")
print("")
print("This script is used for submitting many rendering jobs in parallel. NOTE this submits render jobs for EVERY .pdb trajectory in the current directory")
print('--------------------------------------------------------------------')

def get_ResIDs(traj_name):
    
    drug_name = traj_name.split("_")[7]
    excip_name = traj_name.split("_")[8]
    
    assert os.path.exists(drug_name+".mol2"), "Error! Drug mol2 file missing!"
    assert os.path.exists(excip_name+".mol2"), "Error! Excip mol2 file missing!"
    
    mol2_files = [drug_name+".mol2", excip_name+".mol2"]
    
    resNames = [md.load(filename).top.residue(0).name for filename in mol2_files]
    
    print(drug_name)
    print(excip_name)
    return resNames

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

# Get a list of all .pdb files in the current directory
traj_files = glob.glob('openmm*')

# Iterate through each file
for file in traj_files:
    # Get the systemTag by using the file name
    tag = file.split('.')[0]
    print(f'Writing Files for {tag}\n')
    
    # Get the jobname by splitting the file name on '_' and taking the first 3 letters of the 3rd and 4th fields
    jobname = "rndr_"+ tag.split("_")[-5][:3] +tag.split("_")[-4][:3]+'_'+tag.split("_")[-1]
    res_IDs = get_ResIDs(file)

    writeJob(jobName = jobname, command =f'python render_jawBreaker.py -i {file} -n movies --dr {res_IDs[0]} --er {res_IDs[1]} -f 15 -d')
    print("\n")
print("Job's Done!")

