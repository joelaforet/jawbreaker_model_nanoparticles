import pandas as pd
import numpy as np
import os
from distutils.spawn import find_executable
#from pymol import cmd,stored
import logging
import subprocess
import sys
import argparse
import random
import string
#####################################################################

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=str, help="Input .tsv file that contains name and SMILES of molecules", required=True)
parser.add_argument("-c", "--charge", type = int, default = 0, help = "Input net charge of the input molecule. (Default = 0)")
parser.add_argument("--pH", type=float, default = 7.0, help="Input pH to protonate molecule at.")
args = parser.parse_args()



#####################################################################

def make_mol2(ligand, pH):
    """
    Generates .mol2 and .frcmod files for ligand.

    @param ligand: input ligand name

    Returns
    -------
    .mol2/.frcmod/.pdb files for input ligand.
    """

    smiles = pd.read_csv('{}.tsv'.format(ligand), delimiter = '\t', header = None)[1][0]
    name = ligand


    try:
        print("Running Obabel...")
        os.system(f'obabel -:"{smiles} {name}" -omol2 -p {pH} -r --conformer --gen3d --weighted -h --ff GAFF --partialcharge eem -O {name}.mol2')
        print("Obabel Success!")
    except:
        print("Error! Mol2 file not generated. Check your input files.")
        return

    try:
        print("Running Antechamber...")
        resName = 'Z'.join(random.choice(string.ascii_uppercase) for _ in range(2))
        os.system(f"antechamber -i {name}.mol2 -fi mol2 -o {name}.mol2 -rn {resName} -fo mol2 -c bcc -s 1 -nc {charge}")
        print("Antechamber run successful!")
        os.system(f"obabel {name}.mol2 -O {name}.pdb")
        print("PDB file generated!")
    except:
        print("Antechamber failed, odd number of electrons detected.")
        #return
        print("Running Obabel...")
        os.system("obabel -ipdb {}_h.pdb -omol2 -O {}.mol2".format(ligand,ligand))
        print("Obabel run successful!")

    print("Generating .frcmod file...")
    print("Running parmchk2...")
    try:
        os.system(f"parmchk2 -i {name}.mol2 -f mol2 -o {name}.frcmod -s gaff2 -a Y")
        print("Parmchk2 run successful!")
        print("Ligand successfully parametrized!")
    except:
        print("Error! .frcmod not generated!")

lig = str(args.infile)
pH = args.pH
charge = args.charge

print("Beginning Parametrization...")
ligand = lig.split(".")[0]

print(ligand)

make_mol2(ligand, pH)

print("Job's Done!")



