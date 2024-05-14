'''
make_snowGlobe.py
by Joe Laforet Jr.
jrl78@duke.edu

Following classical nucleation theory, this script generates a spherical nanoparticle of N drug molecules in the center of a box
containing a user specified number (M) of excipient molecules.
If specified, will also solvate system in either PBS or water and generate input files for simulation with OpenMM.

usage: make_snowGlobe.py [-h] -d DRUG -e EXCIP --dNum DNUM --eNum ENUM [-b BUFFER] [--solvent_type {1,2,3,4}] [-# SEED] [--noExcip]
	optional arguments:
	  -h, --help            show this help message and exit
	  -d DRUG, --drug DRUG
				Input name of Drug
      -e EXCIP, --excipient EXCIP
				Input name of Excipient
      --dNum DNUM,
                Input number of Drug molecules to add
      --eNum ENUM,
                Input number of Excipient molecules to add
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

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import defaultdict

import mdtraj as md

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
    #UNITS == NANOMETERS
    name = name
    
    diameter = ( (6*V_molecule)/np.pi)**(1/3)
    print("Diameter of {}: {:.2f} nm".format(name, diameter))

    return diameter

def calc_core_dim(num_drug, v_drug):
    # Solve for radius of a sphere (drug core) given an input volume
    core_volume = num_drug * v_drug #cubic nanometers
    core_radius = ( (3/ (4*np.pi))*core_volume)**(1/3)
    
    print(f"Nanoparticle Diameter: {2*core_radius}nm")
    return np.round(core_radius,2)

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

SHELL_TEMPLATE = """
structure %s.pdb
  number %d
  outside sphere 0. 0. 0. %f
  inside sphere 0. 0. 0. %f
end structure
"""

EXCIP_TEMPLATE = """
structure %s.pdb
  number %d
  inside box %f %f %f %f %f %f
  outside box %f %f %f %f %f %f
end structure
"""

DMSO_TEMPLATE = """
structure %s.pdb
  number %d
  inside box %f %f %f %f %f %f
  outside sphere 0. 0. 0. %f
end structure
"""


def pack_snowGlobe(names, numbers, vols, DMSO = False, box_buffer = 10, tolerance=2.0, seed=-1):
    """Run packmol to generate a snowGlobe containing a mixture of molecules.
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

    for i in range(len(names)):	
        output_filename = output_filename + names[i]+'_'
    output_filename = ('snowGlobe_'+ output_filename + str(numbers[0])+"_" + str(numbers[1])+ '.pdb') 

    rad_core = calc_core_dim(num_drug, vols[0])
    excip_vol = vols[1]
    #need to convert from nm to angstroms
    rad_core*= 10
    excip_dia = dia_from_volume(excip_vol) * 10
    
    drug_section = DRUG_TEMPLATE % (names[0], num_drug, rad_core)

    eBil = (20 + 2*rad_core)/2 # Excipient box inside length
    eBol = eBil + 2*excip_dia # Excipient box outside length

    excip_section = EXCIP_TEMPLATE % (names[1], numbers[1], -eBol, -eBol, -eBol, eBol, eBol, eBol, -eBil, -eBil, -eBil, eBil, eBil, eBil)
    
    snowGlobe_diameter = 2*(rad_core) *0.1 # diameter in nm
    system_volume = (snowGlobe_diameter + 2+ 2*excip_diameter + 2*(buffer_distance/10)) **3 # in nanometers

    header = HEADER_TEMPLATE % (tolerance, output_filename, seed)

    header += drug_section
    if(not noExcip):
        header += excip_section

    if(DMSO):
        num_DMSO = get_DMSO(system_volume)
        print('Num DMSO is: {}'.format(num_DMSO))
        #inside box %f %f %f %f %f %f
        #outside sphere 0. 0. 0. %f

        DMSO_section = DMSO_TEMPLATE % ('DMSO', num_DMSO, -eBil, -eBil, -eBil, eBil, eBil, eBil, rad_core + tolerance + 12)
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
            return output_filename, snowGlobe_diameter
        else:
            os.system("%s < %s" % (PACKMOL_PATH, file_handle.name))
            return output_filename, snowGlobe_diameter
            # Returns files name of jawBreaker and diameter of full system in Angstroms
    except:
        print("PackMol Failed!")


# MAIN CODE	
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-d", "--drug", type=str, required=True, help="Input name of Drug")
parser.add_argument("-e", "--excip", type=str, required=True, help="Input name of Excipient")
parser.add_argument("--dNum", type=int, required=True, help="Input number of drug molecules to add")
parser.add_argument("--eNum", type=int, required=True, help="Input number of excipient molecules to add")
parser.add_argument('--solvent_type', type=int, default=1, choices=range(1, 6), help="1 = Vacuum ; 2 = TIP3P ; 3 = PBS; 4 = PBS+DMSO ; 5 = Implicit PBS+DMSO")
parser.add_argument("-b", "--buffer", type=int, default=10, help="Buffering distance between waterbox and molecule box (default=10 Angstroms)")
parser.add_argument("-#", "--seed", type=int, default=-1, help="If specified, creates a configuration with given seed number (default=random)")
parser.add_argument("--noExcip", type = bool, default = False, help="If specified, generates a snowGlobe with no excipient shell.")
args = parser.parse_args()

print('--------------------------------------------------------------------')
print("     Welcome to make_snowGlobe.py!      ")
print("     By Joe Laforet Jr.      ")
print("")
print("This script is used for constructing Drug:Excipient Nanoparticles following the SnowGlobe hypothesis.")
print('--------------------------------------------------------------------')

names = [args.drug, args.excip]
numbers = [args.dNum, args.eNum]
resNames = []

num_drug = numbers[0]
num_excip = numbers[1]
DMSO = False

noExcip = args.noExcip

if args.solvent_type == 4 or args.solvent_type == 5:
    DMSO = True

print("Running initial calculations...")

drug_smiles = pd.read_csv('{}.tsv'.format(names[0]), delimiter = '\t', header = None)[1][0]
excip_smiles = pd.read_csv('{}.tsv'.format(names[1]), delimiter = '\t', header = None)[1][0]

v_drug = smiles_to_vol(drug_smiles)
v_excip = smiles_to_vol(excip_smiles)

excip_diameter = dia_from_volume(v_excip)

vols = [v_drug, v_excip]

rad_core = calc_core_dim(num_drug, vols[0])

excip_box_length = 2*rad_core + 20 # in angstroms


if(noExcip):
    num_excips = 0

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
addions box Cl- 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box PBS_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

TLEAP_PBS_DMSO_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
STK = loadmol2 DMSO.mol2
box = loadPdb %(box_filename)s
%(amberparams_section)s
loadamberparams DMSO.frcmod
addions box Na+ 0
addions box Cl- 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box PBS-DMSO_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

TLEAP_IMP_PBS_DMSO_TEMPLATE = """
source leaprc.gaff2
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
STK = loadmol2 DMSO.mol2
box = loadPdb %(box_filename)s
%(amberparams_section)s
loadamberparams DMSO.frcmod
addions box Na+ 0
addions box Cl- 0
addions box Na+ %(num_Sodium)s
addions box K+ %(num_Potassium)s
addions box Cl- %(num_Chlorine)s
setbox box centers
savepdb box IMP_PBS-DMSO_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""

resNames = [md.load(filename).top.residue(0).name for filename in mol2_files]

print('--------------------------------------------------------------------')
print("Making the SnowGlobe!")

outName, snowGlobe_diameter = pack_snowGlobe(names, numbers, vols, DMSO = DMSO, box_buffer = buffer_distance, seed=seed)
#diameter in nanometers
if(exists(outName)):
    print("SnowGlobe made successfully.")
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

    system_volume = (snowGlobe_diameter*10 + 20 + 2*(excip_diameter*10) + 2*(buffer_distance)) **3 # cubic angstroms

    print('SnowGlobe diameter is {}'.format(snowGlobe_diameter))
    print('Buffer distance is {}'.format(buffer_distance))
    print('System volume is {}'.format(system_volume))
    #For these systems we are only going to consider the effects of NaCl and KCl
    #The polyatomic ions are present in negligible amounts at our current simulation size,
    #but in the future they can also be incorporated.
    for x in range(2):
        ions_needed.append(round(get_molecules(system.iloc[x][1], system_volume)))

    print(ions_needed)
    tleap_commands = TLEAP_PBS_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, num_Sodium = ions_needed[0],
                                             num_Potassium = ions_needed[1], num_Chlorine = ions_needed[0] + ions_needed[1], buffer_distance=buffer_distance,
                                            box_filename=outName, prmtop_filename='PBS_'+outName.split('.')[0] + '.prmtop', inpcrd_filename='PBS_' +outName.split('.')[0] +'.inpcrd')

elif args.solvent_type == 4: #TIP3P water with PBS+1% DMSO:

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

    system_volume_nm = (snowGlobe_diameter + 2+ 2*excip_diameter + 2*(buffer_distance/10)) **3 # in nm

    system_volume_ang = (snowGlobe_diameter*10 + 20 + 20*excip_diameter + 2*(buffer_distance)) **3 # in ang

    print('SnowGlobe diameter is {}'.format(snowGlobe_diameter))
    print('Buffer distance is {}'.format(buffer_distance))
    print('System volume is {}'.format(system_volume_nm))

    #For these systems we are only going to consider the effects of NaCl and KCl
    #The polyatomic ions are present in negligible amounts at our current simulation size,
    #but in the future they can also be incorporated.
    for x in range(2):
        ions_needed.append(round(get_molecules(system.iloc[x][1], system_volume_ang)))

    print(ions_needed)
    
    #print(DMSO_mols)
    tleap_commands = TLEAP_PBS_DMSO_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, num_Sodium = ions_needed[0], num_Potassium = ions_needed[1], num_Chlorine = ions_needed[0] + ions_needed[1], buffer_distance=buffer_distance, box_filename=outName, prmtop_filename='PBS-DMSO_'+outName.split('.')[0] + '.prmtop',inpcrd_filename='PBS-DMSO_' +outName.split('.')[0] +'.inpcrd')

elif args.solvent_type == 5: #Implicit water with PBS+1% DMSO:

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

    system_volume_nm = (snowGlobe_diameter + 2 + 2*excip_diameter + 2*(buffer_distance/10)) **3 # in nm

    system_volume_ang = (snowGlobe_diameter*10 + 20 + 20*(excip_diameter) + 2*(buffer_distance)) **3 # cubic angstroms


    print('SnowGlobe diameter is {}'.format(snowGlobe_diameter))
    print('Buffer distance is {}'.format(buffer_distance))
    print('System volume is {}'.format(system_volume_nm))

    #For these systems we are only going to consider the effects of NaCl and KCl
    #The polyatomic ions are present in negligible amounts at our current simulation size,
    #but in the future they can also be incorporated.
    for x in range(2):
        ions_needed.append(round(get_molecules(system.iloc[x][1], system_volume_ang)))

    print(ions_needed)
    
    #print(DMSO_mols)
    tleap_commands = TLEAP_IMP_PBS_DMSO_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, num_Sodium = ions_needed[0], num_Potassium = ions_needed[1], num_Chlorine = ions_needed[0] + ions_needed[1], box_filename=outName, prmtop_filename='IMP_PBS-DMSO_'+outName.split('.')[0] + '.prmtop',inpcrd_filename='IMP_PBS-DMSO_' +outName.split('.')[0] +'.inpcrd')





file_handle = open('tleap_commands', 'w')
file_handle.writelines(tleap_commands)
file_handle.close()
print('Success!')
print('--------------------------------------------------------------------')
logger.debug('Running tleap in temporary directory.')
cmd = "tleap -f %s " % file_handle.name
logger.debug(cmd)

output = getoutput(cmd)
logger.debug(output)


print("Job's Done!")
print('--------------------------------------------------------------------')

