'''
make_jawBreaker.py
by Joe Laforet Jr.
jrl78@duke.edu

Generates spherical nanoparticle consisting of N drug molecule core surrounded by a shell of excipient.
If specified, will also solvate system in either PBS or water and generate input files for simulation with OpenMM.

usage: make_jawBreaker.py [-h] -d DRUG -e EXCIP --dNum DNUM [-b BUFFER] [--solvent_type {1,2,3,4}] [-# SEED] [--noExcip]
	optional arguments:
	  -h, --help            show this help message and exit
	  -d DRUG, --drug DRUG
				Input name of Drug
      -e EXCIP, --excipient EXCIP
				Input name of Excipient
      --dNum DNUM,
                Input number of Drug molecules to add
	  -b BUFFER, --buffer BUFFER
				Buffering distance between waterbox and molecule box
				(default=10 Angstroms)
      --solvent_type {1,2,3,4,5}
                1 = Vaccuum (default); 2 = TIP3P Water ; 3 = PBS; 4 = PBS+DMSO; 5 = Implicit PBS+DMSO
      -# SEED, --seed SEED  If specified, creates a configuration with given seed
				number (default=random)
      --noExcip,
                If specified, generates a jawBreaker with no excipient molecules. (default = False)
            
'''

import pandas as pd
import numpy as np

import argparse

import os
from os.path import exists
from distutils.spawn import find_executable

import math

import logging
import subprocess
import sys

import re

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import defaultdict

import mdtraj as md
import requests

logger = logging.getLogger(__name__)


def getoutput(cmd):
    """Compatibility function to substitute deprecated commands.getoutput in Python2.7 (Original code : openmoltools.amber)""" 
    try:
        out = subprocess.getoutput(cmd)
    except AttributeError:
        out = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT,
                               stdout=subprocess.PIPE).communicate()[0]
    try:
        return str(out.decode())
    except:
        return str(out)

def smiles_to_vol(smiles):
    """Calculate the volume of a molecule given its SMILES using the VdW approximation.
    Parameters:
    ----------------------------------------
    smiles: (str) SMILES from online resource
    
    Returns: (float) Volume of the molecule in cubic Nm
    """
    
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    formula = CalcMolFormula(mol)

    atomic_count = defaultdict(lambda : 0)
    
    for atom in mol.GetAtoms():
        atomic_count[atom.GetSymbol()] += 1
            
    total_vol = 0
    
    #Dictionary contains Atomic Symbol and VdW radius retrieved from Wikipedia
    atoms_to_parse = {'C':1.7 , 'H':1.2, 'N':1.55, 'O':1.52, 'P':1.8, 'Cl':1.75, 'S':1.8, 'F':1.47, 'I':1.98}
    
    for x in atoms_to_parse.keys():
        if(x not in atoms_to_parse.keys()):
            raise Exception("Error, Atom not found in VdW dictionary. Find VdW for {} and add it to dictionary.".format(x))
        
        total_vol += atomic_count[x] * (4*np.pi* (atoms_to_parse[x]**3))/3
    
    total_vol/=1000
    #print("Volume of molecule is: {:.4f} Cubic Nm".format(total_vol))
    
    return total_vol

def dia_from_volume(V_molecule, name = 'Molecule'):
    #We calculate the "diameter" of one molecule based on its VdW volume again
    #This is used to make the single diameter thick excipient coating.

    diameter = ( (6*V_molecule)/np.pi)**(1/3)
    print("Diameter of {}: {:.2f} nm".format(name, diameter))

    return diameter

def calc_core_dim(num_drug, v_drug):
    # Solve for radius of a sphere (drug core) given an input volume
    core_volume = num_drug * v_drug
    core_radius = ( (3/ (4*np.pi))*core_volume)**(1/3)
    
    return core_volume, core_radius

def calc_sphere_dims(num_drug, v_drug, v_excip, names = ['Drug', 'Excip']):
    #Calculate the dimensions of the jawBreaker. Returns radius of core, diameter of one excipient, and the number of excipients needed to fill the shell.

    excip_dia = dia_from_volume(v_excip, names[1])
    
    vol_core, rad_core = calc_core_dim(num_drug, v_drug)
    
    vol_shell = (4/3)*np.pi * ( (rad_core + excip_dia)**3 - rad_core**3)
    
    num_excips = vol_shell / v_excip
    
    num_excips += 0 * num_excips # implementing the fudge factor
    if(noExcip):
        num_excips = 0

    print("Nanoparticle has Core Radius: {:.2f}nm\nExcipients needed to fill a shell: {:.2f}.".format(rad_core, num_excips))
    
    return np.round(rad_core,2) , np.round(excip_dia,2), np.round(num_excips)

def ang_to_L(ang):
    return ang*1e-27

def get_molecules(conc, vol):
    
    #Concentration in Molar (mols/L)
    #Volume provided in cubic angstrom
    
    new_vol = ang_to_L(vol)
    
    return conc *new_vol*6.02e23

def get_DMSO(vol):

    #We assume the solution is prepared with 1% DMSO in volume ratio
    fraction = vol * 0.01
    DMSO_vol = 0.1237 # nm^3 by VdW approximation

    return round(fraction/DMSO_vol)



PACKMOL_PATH = find_executable("packmol")

HEADER_TEMPLATE = """
# Mixture 
tolerance %f
filetype pdb
output %s
add_amber_ter
seed %d
"""

DRUG_TEMPLATE = """
structure %s.pdb
  number %d
  inside sphere 0. 0. 0. %f
end structure
"""

EXCIP_TEMPLATE = """
structure %s.pdb
  number %d
  outside sphere 0. 0. 0. %f
  inside sphere 0. 0. 0. %f
end structure
"""

def pack_Sphere(names, numbers, DMSO = False, box_buffer = 10, tolerance=2.0, seed=-1):
    """Run packmol to generate a sphere containing a mixture of molecules.
    Parameters
    ----------
    names : list(str)
        List of pdb filenames. Omit the file extension. i.e. [Sorafenib, Glycyrrhizin] 
    numbers : list(int)
        The number of molecules of each mixture component. i.e. [#Drug, #Excip]
    DMSO : bool (default: False)
        If specified, will add 1% DMSO molecules to box
    box_buffer : float (default = 10.0)
        Buffer distance of water box around jawBreaker. (Angstroms)
    tolerance : float, optional, default=2.0
        The minimum spacing between molecules during packing.  In ANGSTROMS!
    seed : default = -1 (random). User may specify seed number for consistent configuration generation
    ----- 
    """

    if PACKMOL_PATH is None:
        raise(IOError("Packmol not found, cannot run pack_Sphere()"))

    # The path to packmol's output PDB file. Concatenates name and number of molecules
    if(DMSO):
        output_filename = 'DMSO_'
    else:
        output_filename = ''

    #for i in range(len(names)):	
     #   output_filename = output_filename + names[i]+'_'

    output_filename = output_filename + names[0] + '_' + names[1] + '_'
    output_filename = (f'{buffer_distance}A_jawBreaker_'+ output_filename + str(numbers[0])+"_" + str(numbers[1])+ '.pdb') 

    rad_core, excip_dia, num_excips = calc_sphere_dims(num_drug, v_drug, v_excip, names)
    
    #need to convert from nm to angstroms
    rad_core*= 10
    excip_dia *=10
    
    drug_section = DRUG_TEMPLATE % (names[0], num_drug, rad_core)
    excip_section = EXCIP_TEMPLATE % (names[1], num_excips, rad_core, rad_core+excip_dia)
    
    jawBreaker_diameter = 2*(rad_core + excip_dia) *0.1 # diameter in nm
    #print('JawBreaker Diameter is: {}'.format(jawBreaker_diameter))

    system_volume = (jawBreaker_diameter + 2*(box_buffer/10)) **3 # in nanometers

    #print('System volume is: {}'.format(system_volume)) # volume in nm^3 ?

    header = HEADER_TEMPLATE % (tolerance, output_filename, seed)

    header += drug_section
    if(not noExcip):
        header += excip_section

    if(DMSO):
        num_DMSO = get_DMSO(system_volume)
        print('Num DMSO is: {}'.format(num_DMSO))
        DMSO_section = EXCIP_TEMPLATE % (names[2], num_DMSO, rad_core+excip_dia, rad_core+excip_dia + tolerance + 12)
        header += DMSO_section


	# Create input file for packmol.
    packmol_filename = output_filename+'.inp'
    with open(packmol_filename, 'w') as file_handle:
        file_handle.write(header)
    print('--------------------------------------------------------------------')
    print("Running PackMol...")
	# Run packmol.

    try:
        if exists(output_filename):
            print("System file detected! Proceeding with Parametrization.")
            return output_filename, jawBreaker_diameter
        else:
            os.system("%s < %s" % (PACKMOL_PATH, file_handle.name))
            return output_filename, jawBreaker_diameter
            # Returns files name of jawBreaker and diameter of full system in Angstroms
    except:
        print("PackMol Failed!")


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

# Function to calculate the amount of ions needed for PBS
def get_Ions(jawBreaker_diameter, buffer_distance):
    #Salt Concentrations for PBS in Molar (mols/L)

    """
    NaCl = .137
    KCl = 0.0027
    Na2HPO4 = 0.008
    KH2PO4 = 0.002
    """

    system = pd.DataFrame()
    system['Ion'] = ['NaCl', "KCl", "Na2HPO4", "KH2PO4"]

    system['Concentration (M)'] = [0.137, 0.0027, 0.008, 0.002]

    ions_needed = []

    system_volume = (jawBreaker_diameter*10 + 2*(buffer_distance)) **3 # cubic angstroms

    print('JawBreaker diameter is {}'.format(jawBreaker_diameter))
    print('Buffer distance is {}'.format(buffer_distance))
    print('System volume is {}'.format(system_volume))

    #For these systems we are only going to consider the effects of NaCl and KCl
    #The polyatomic ions are present in negligible amounts at our current simulation size,
    #but in the future they can also be incorporated.

    for x in range(2):
        ions_needed.append(round(get_molecules(system.iloc[x][1], system_volume)))

    print(ions_needed)
    return ions_needed


# MAIN CODE	
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-d", "--drug", type=str, required=True, help="Input name of Drug")
parser.add_argument("-e", "--excip", type=str, required=True, help="Input name of Excipient")
parser.add_argument("--DMSO", type=str, required=True, help="Input name of DMSO (for alternate charge schemes)")
parser.add_argument("--dNum", type=int, required=True, help="Input number of drug molecules to add")
parser.add_argument('--solvent_type', type=int, default=1, choices=range(1, 6), help="1 = Vacuum ; 2 = TIP3P ; 3 = PBS; 4 = PBS+DMSO ; 5 = Implicit PBS+DMSO")
parser.add_argument("-b", "--buffer", type=int, default=10, help="Buffering distance between waterbox and molecule box (default=10 Angstroms)")
parser.add_argument("-#", "--seed", type=int, default=-1, help="If specified, creates a configuration with given seed number (default=random)")
parser.add_argument("--noExcip", type = bool, default = False, help="If specified, generates a jawBreaker with no excipient shell.")
args = parser.parse_args()

print('--------------------------------------------------------------------')
print("     Welcome to make_jawBreaker.py!      ")
print("     By Joe Laforet Jr.      ")
print("")
print("This script is used for constructing Drug:Excipient Nanoparticles following the Jawbreaker hypothesis.")
print('--------------------------------------------------------------------')

names = [args.drug, args.excip, args.DMSO]
numbers = [args.dNum]
resNames = []

num_drug = numbers[0]

DMSO = False

noExcip = args.noExcip

if args.solvent_type == 4 or args.solvent_type == 5:
    DMSO = True

print("Running initial calculations...")

#Function to ensure that query inputs are in pubchem parsable format
# I.E. CholicAcid -> Cholic Acid

split_string = lambda s: re.sub(r'(?<=\w)(?=[A-Z])', ' ', s)


#drug_smiles = pd.read_csv('{}.tsv'.format(names[0]), delimiter = '\t', header = None)[1][0]
drug_smiles = get_SMILES(split_string(names[0]))
excip_smiles = get_SMILES(split_string(names[1]))
#excip_smiles = pd.read_csv('{}.tsv'.format(names[1]), delimiter = '\t', header = None)[1][0]

v_drug = smiles_to_vol(drug_smiles)
v_excip = smiles_to_vol(excip_smiles)

rad_core, excip_dia, num_excips = calc_sphere_dims(num_drug, v_drug, v_excip, names)
if(noExcip):
    num_excips = 0

numbers.append(round(num_excips))

buffer_distance = args.buffer # angstroms

seed=args.seed
print('--------------------------------------------------------------------')
# Check for one residue name per mol2 file and uniqueness between all mol2 files
all_names = set()

mol2_files = [x+'.mol2' for x in names]
frcmod_files = [x + '.frcmod' for x in names]

for filename in mol2_files:
    t = md.load(filename)
    inFile_names = set([r.name for r in t.top.residues])

    if len(inFile_names) != 1:
        raise(ValueError("Must have a SINGLE residue name in each mol2 file."))

    resNames.append(list(inFile_names)[0])
    all_names = all_names.union(list(inFile_names))

if len(all_names) != len(mol2_files):
    raise(ValueError("Must have UNIQUE residue names in each mol2 file."))
if len(mol2_files) != len(frcmod_files):
    raise(ValueError("Must provide an equal number of frcmod and mol2 file names."))

# Checks passed, now to fill the parametrization template

TLEAP_VAC_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Cl- 0
setbox box centers
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

TLEAP_TIP3P_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Cl- 0
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box WAT_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

TLEAP_PBS_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
addions box Cl- 0
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box PBS_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

#Note for this run I removed the line:
#STK = loadmol2 DMSO.mol2
TLEAP_PBS_DMSO_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
loadamberparams DMSO.frcmod
addions box Na+ 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
addions box Cl- 0
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box PBS-DMSO_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

#Note for this run I removed the line:
#STK = loadmol2 DMSO.mol2
TLEAP_IMP_PBS_DMSO_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
loadamberparams DMSO.frcmod
addions box Na+ 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
addions box Cl- 0
setbox box centers
savepdb box IMP_PBS-DMSO_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

resNames = [md.load(filename).top.residue(0).name for filename in mol2_files]

print('--------------------------------------------------------------------')
print("Making the JawBreaker!")

outName, jawBreaker_diameter = pack_Sphere(names, numbers, DMSO = DMSO, box_buffer = buffer_distance, seed=seed)
#diameter in nanometers
if(exists(outName)):
    print("JawBreaker made successfully.")
else:
    print("Uh Oh! Something went wrong.")
    print("See above error from PackMol!")
    exit()

print('--------------------------------------------------------------------')
print("Solvating and Parametrizing with tleap...")

mol2_section = "\n".join("%s = loadmol2 %s" % (
        resNames[k], filename) for k, filename in enumerate(mol2_files))

amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(frcmod_files))

if args.solvent_type == 1: #Vacuum
    tleap_commands = tleap_commands = TLEAP_VAC_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section,
                                            box_filename=outName, prmtop_filename='VAC_'+outName.split('.')[0] + '.prmtop', inpcrd_filename='VAC_' +outName.split('.')[0] +'.inpcrd')
elif args.solvent_type == 2: #Just TIP3P water:
    tleap_commands = TLEAP_TIP3P_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, buffer_distance=buffer_distance,
                                            box_filename=outName, prmtop_filename='WAT_'+outName.split('.')[0] + '.prmtop', inpcrd_filename='WAT_' +outName.split('.')[0] +'.inpcrd')
elif args.solvent_type == 3: #TIP3P water with PBS:
    ions_needed = get_Ions(jawBreaker_diameter = jawBreaker_diameter, buffer_distance=buffer_distance)

    tleap_commands = TLEAP_PBS_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, num_Sodium = ions_needed[0],
                                             num_Potassium = ions_needed[1], num_Chlorine = ions_needed[0] + ions_needed[1], buffer_distance=buffer_distance,
                                            box_filename=outName, prmtop_filename='PBS_'+outName.split('.')[0] + '.prmtop', inpcrd_filename='PBS_' +outName.split('.')[0] +'.inpcrd')

elif args.solvent_type == 4: #TIP3P water with PBS+1% DMSO:
    ions_needed = get_Ions(jawBreaker_diameter = jawBreaker_diameter, buffer_distance=buffer_distance)

    tleap_commands = TLEAP_PBS_DMSO_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section,
                                                    num_Sodium = ions_needed[0], num_Potassium = ions_needed[1],
                                                    num_Chlorine = ions_needed[0] + ions_needed[1], buffer_distance=buffer_distance,
                                                    box_filename=outName, prmtop_filename='PBS-DMSO_'+outName.split('.')[0] + '.prmtop',
                                                    inpcrd_filename='PBS-DMSO_' +outName.split('.')[0] +'.inpcrd')

elif args.solvent_type == 5: #Implicit water with PBS+1% DMSO:
    ions_needed = get_Ions(jawBreaker_diameter = jawBreaker_diameter, buffer_distance=buffer_distance)

    tleap_commands = TLEAP_IMP_PBS_DMSO_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section,
                                                        num_Sodium = ions_needed[0], num_Potassium = ions_needed[1],
                                                        num_Chlorine = ions_needed[0] + ions_needed[1], box_filename=outName,
                                                        prmtop_filename='IMP_PBS-DMSO_'+outName.split('.')[0] + '.prmtop',
                                                        inpcrd_filename='IMP_PBS-DMSO_' +outName.split('.')[0] +'.inpcrd')


file_handle = open('tleap_commands', 'w')
file_handle.writelines(tleap_commands)
file_handle.close()
print('Success!')
print('--------------------------------------------------------------------')
logger.debug('Running tleap in temporary directory.')
cmd = f"tleap -f {file_handle.name} -o {outName.split('.')[0]}.log"
logger.debug(cmd)

output = getoutput(cmd)
logger.debug(output)


print("Job's Done!")
print('--------------------------------------------------------------------')
