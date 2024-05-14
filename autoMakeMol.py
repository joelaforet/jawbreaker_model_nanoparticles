import glob
import subprocess

files = glob.glob("*.tsv")
for file in files:
    molName = file.split('.')[0]
    command = f'python slurm.py -i "python make_mol2.py -i {molName} --pH 7.4" -r -n "make{molName}"'
    subprocess.run(command, shell=True)
print("Job's Done!")


