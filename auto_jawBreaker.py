import glob
import subprocess

# Get a list of all .mol2 files in the current directory
files = glob.glob('*.mol2')

# Iterate through each file
for file in files:
    # Split the file name on '.' and take the first field as the excipient
    excipient = file.split('.')[0]
    
    #exclude_list = ["Indomethacin", "CongoRed", "EvansBlue", "Glycyrrhizin", "Ursodiol", "Budesonide", "Candesartan", "FolicAcid", "Sorafenib", "Carbamazepine"]
    
    #if excipient in exclude_list:
     #   continue
    # Construct the command with the excipient


    command = f'python slurm.py -i "python make_jawBreaker.py -d Fulvestrant -e {excipient} --dNum 25 --solvent_type 3 -b 15" -r -n "jb{excipient}"'    
   # Execute the command using the subprocess module
    subprocess.run(command, shell=True)
