import os
import subprocess
import pandas as pd

import argparse
import random
import string
import requests
import shutil

from openbabel import openbabel as ob

###################################################################

# Define the get_SMILES function
def get_SMILES(query):
    # Construct the PubChem API URL with the search query
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/property/CanonicalSMILES/JSON"

    # Send a GET request to the API URL
    response = requests.get(url)

    # Check if the response was successful
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()

        # Extract the canonical SMILES string from the JSON response
        smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]

        # Print the canonical SMILES string
        print(f"Canonical SMILES string for {query}: {smiles}")
        return smiles
    else:
        # Print an error message if the response was not successful
        print(f"Error: {response.status_code}")

def calculate_net_charge(smiles):
    # Create an OpenBabel molecule object from the SMILES string
    mol = ob.OBMol()
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("smiles")
    obConversion.ReadString(mol, smiles)

    # Calculate the total charge of the molecule
    total_charge = 0
    for atom in ob.OBMolAtomIter(mol):
        total_charge += atom.GetFormalCharge()

    return total_charge


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
        fh.write("#SBATCH --mem=20G\n")
        fh.write("#SBATCH -p scavenger-gpu --gres=gpu:1\n")
        fh.write("#SBATCH --exclusive\n")
        fh.write('%s\n' % command)
    # Submit the job using sbatch
    subprocess.run(['sbatch', job_file])

######################################################################################

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=str, help="Input .csv file that contains name and of molecules", required=True)
#parser.add_argument("-c", "--charge", type = int, default = 0, help = "Input net charge of the input molecule. (Default = 0)")
args = parser.parse_args()

#################################################################################
to_make = args.infile


df = pd.read_csv(to_make, header = None).squeeze("columns")

mols = df.tolist()

# Define the list of drugs to convert
#drugs = ['aspirin', 'ibuprofen', 'acetaminophen']

# Loop through each drug and submit a job to convert to mol2
for molecule in mols:
    # Get the SMILES string for the drug
    smiles = get_SMILES(molecule)
    molecule = "".join(molecule.split())
    charge = calculate_net_charge(smiles)

    command = 'INPUT_DIR="."\nOUTPUT_DIR="."\n'

    command+= f'qmmm_dir={molecule}_qmmm\n'
    command+= f'mkdir "$qmmm_dir"\n'
    # Command to convert SMILES to mol2
    command += f'obabel -:"{smiles} {molecule}" -omol2 -r --conformer --gen3d --weighted -h --partialcharge eem -O {molecule}.mol2\n'
    command+= f'mv {molecule}.mol2 {molecule}_qmmm\ncd {molecule}_qmmm\n'

    # Add a unique residue ID to the mol2 and calculate partial charges w/ Antechamber
    resName = 'Z' + ''.join(random.choice(string.ascii_uppercase) for _ in range(2))

    command += "\n" + f"antechamber -i {molecule}.mol2 -fi mol2 -o {molecule}.mol2 -rn {resName} -fo mol2 -s 1 -nc {charge}"


    # Write a pdb file for the mol2 file
    command += "\n" + f"antechamber -i {molecule}.mol2 -fi mol2 -o {molecule}.pdb -fo pdb -dr y"

    # Write a frcmod file for the mol2 file
    command += "\n" + f"parmchk2 -i {molecule}.mol2 -f mol2 -o {molecule}.frcmod -s gaff2 -a Y"

    command += "\n" +f"""# Move the output files to the original directory
    for output_file in *; do
        if [[ $output_file != "sqm.in" && $output_file != "sqm.out" && $output_file != "sqm.pdb" ]]; then
            mv "$output_file" "../"
        fi
    done

    # Delete the QM/MM directory and its remaining files
    cd "../"\n
    rm -r {molecule}_qmmm\n"""


    # Write and submit the job to the cluster
    writeJob(f"make_{molecule}", command)

print("Job's Done!")
