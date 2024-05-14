"""
manySims.py
By Joe Laforet Jr.
jrl78@duke.edu

This script serves as an easy way to queue multiple simulations at once. A user can specify the amount of replicates to run, time to simulate, and amount of frames to save. The script will then make, submit, and execute the simulation files in parallel.

usage: many_Sims.py [-h] -r REPLICATES 
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
import argparse
import os

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-r","--runs", type=int, required=True, help="Input number of replicates to run for each simulation")
parser.add_argument("-t","--time", type = int, required = True, help="Input time to simulate molecules for. (Nanoseconds)")
parser.add_argument("-f","--frames", type = int, required = True, help="Input amount of frames to save.")
args = parser.parse_args()

print('--------------------------------------------------------------------')
print("     Welcome to many_Sims.py!      ")
print("     By Joe Laforet Jr.      ")
print("")
print("This script is used for submitting many simulation jobs in parallel, for the desired amount of replicates. NOTE this submits simulations for EVERY .prmtop/.inpcrd pair in the current directory")
print('--------------------------------------------------------------------')



# Get a list of all .pdb files in the current directory
prmtops = glob.glob('*.prmtop')
inpcrds = glob.glob("*.inpcrd")

assert len(prmtops) == len(inpcrds), "Error! Unmatched number of .prmtop and .inpcrd files"

names = [x.split('.')[0] for x in prmtops]
runs = args.runs
frames = args.frames
time = args.time

# Iterate through each file
for file in prmtops:
    # Get the systemTag by using the file name
    tag = file.split('.')[0]
    print(f'Writing Files for {tag}\n')
    for x in range(1,runs+1):
    # Get the jobname by splitting the file name on '_' and taking the first 3 letters of the 3rd and 4th fields
        jobname = tag.split("_")[-4][:3] +tag.split("_")[-2] + tag.split("_")[-3][:3] + tag.split("_")[-1].split(".")[0] + f"r{x}"


        os.system(f"python slurm.py -i 'python run_NPT_openmm_simulation.py -s {tag} --time {time} --frames {frames} -r {x}' -r -n {jobname}")
    print("\n")
print("Job's Done!")

