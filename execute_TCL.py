"""
execute_TCL.py
By Joe Laforet Jr.
jrl78@duke.edu

This code is used for executing all .tcl analysis scripts in a directory at once. This will
run the scripts in series, so beware if the scripts take long. We also assume that vmd
is installed correctly on your current environment.
This code will execute all .tcl scripts in parallel using SLURM on the DCC. *Requires slurm.py to be in same directory.

"""


import glob
import os

# Find all files in the current directory that end in .tcl
tcl_files = glob.glob("*.tcl")

# Loop through each .tcl file and execute the command "vmd < FILE.tcl"
for tcl_file in tcl_files:
    name = tcl_file.split(".")[0]
    os.system(f"python slurm.py -i 'vmd < {tcl_file}' -r -n {name}")

