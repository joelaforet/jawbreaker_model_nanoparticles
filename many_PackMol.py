import os
import glob

files = glob.glob("*.inp")

for file in files:
    os.system(f'python slurm.py -i "packmol < {file}" -r -n {file}')

print("Job's Done")
