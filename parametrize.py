"""

Adapted from make_Box.py
by Joe Laforet
jrl78


"""
'''
Generate AMBER-recognized input files from a tab-delimited file consisting of the molecule components and the number of each molecules per line, as well as a .pdb of a nanoparticle (with antibody) system.
Example line entry of tab-delimited file: 'Sorafenib  CholicAcid    2   2'
usage: make_box.py [-h] -i INFILE [-s] [-b BUFFER] [-# SEED] -n PDB
	optional arguments:
	  -h, --help            show this help message and exit
	  -i INFILE, --infile INFILE
				Input file that contains names and number of molecules
	  -s, --solvate         If activated, creates a model solvated in explicit water (TIP3P)
	  -b BUFFER, --buffer BUFFER
				Buffering distance between waterbox and molecule box
				(default=10 Angstroms)
	  -# SEED, --seed SEED  If specified, creates a configuration with given seed
				number (default=random)
          -n PDB    Input .pdb of antibody and nanoparticle that will be parametrized for OpenMM input.
          -o name   Name for output .prmtop and .inpcrd files.
'''

import mdtraj as md
import parmed
import pandas as pd
import numpy as np
import os
from mdtraj.utils.delay_import import import_
from distutils.spawn import find_executable
import logging
import subprocess
import sys
import argparse

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

def mol2_to_pdb(mol2_filename):
	#convert PDB file with this name, to be used for Packmol input
	struct = md.load_mol2(mol2_filename)
	struct.save_pdb(mol2_filename[:-4]+'pdb')


TLEAP_TEMPLATE = """
source leaprc.gaff
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

TLEAP_TEMPLATE_SOLV = """
source leaprc.gaff
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
%(mol2_section)s
box = loadPdb %(box_filename)s
%(amberparams_section)s
addions box Na+ 0
addions box Cl- 0
solvatebox box TIP3PBOX %(buffer_distance)d
setbox box centers
savepdb box solv_%(box_filename)s
saveAmberParm box %(prmtop_filename)s %(inpcrd_filename)s
quit
"""



def build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, prmtop_filename, inpcrd_filename, solvation=False, buffer_distance=10):
    """Create a prmtop and inpcrd from a collection of mol2 and frcmod files
    as well as a single box PDB.  We have used this for setting up
    simulations of binary mixtures. (Original code : openmoltools.amber)
    Parameters
    ----------
    mol2_filenames : list(str)
        Filenames of GAFF flavored mol2 files.  Each must contain exactly
        ONE ligand.
    frcmod_filenames : str
        Filename of input GAFF frcmod filenames.
    box_filename : str
        Filename of PDB containing an arbitrary box of the mol2 molecules.
    prmtop_filename : str
        output prmtop filename.  Should have suffix .prmtop
    inpcrd_filename : str
        output inpcrd filename.  Should have suffix .inpcrd
    solvation : Boolean, optional. Default: False
        Boolean for whether the system should be solvated explicitly or not. If true, the system will be solvated in TIP3P water model using tleap
    buffer_distance : int, optional. Default: 10
        If solvation is true, will add water molecules with a buffering distance of 10 Angstrom unless specified otherwise.
    Returns
    -------
    tleap_commands : str
        The string of commands piped to tleap for building the prmtop
        and inpcrd files.  This will *already* have been run, but the
        output can be useful for debugging or archival purposes. However,
        this will reflect temporary file names for both input and output
        file as these are used to avoid tleap filename restrictions.
    Notes
    -----
    This can be easily broken if there are missing, duplicated, or
    inconsistent ligand residue names in your box, mol2, and frcmod files.
    You can use mdtraj to edit the residue names with something like
    this: trj.top.residue(0).name = "L1"
    """

    # Check for one residue name per mol2 file and uniqueness between all mol2 files
    all_names = set()
    for filename in mol2_filenames:
        t = md.load(filename)
        names = set([r.name for r in t.top.residues])

        if len(names) != 1:
            raise(ValueError("Must have a SINGLE residue name in each mol2 file."))

        all_names = all_names.union(list(names))

    if len(all_names) != len(mol2_filenames):
        raise(ValueError("Must have UNIQUE residue names in each mol2 file."))
    if len(mol2_filenames) != len(frcmod_filenames):
        raise(ValueError("Must provide an equal number of frcmod and mol2 file names."))

    #Get number of files
    nfiles = len(mol2_filenames)

    #Build absolute paths of output files so we can copy them back
    prmtop_filename = os.path.abspath(prmtop_filename)
    inpcrd_filename = os.path.abspath(inpcrd_filename)

    all_names = [md.load(filename).top.residue(0).name for filename in mol2_filenames]

    mol2_section = "\n".join("%s = loadmol2 %s" % (
        all_names[k], filename) for k, filename in enumerate(mol2_filenames))
    #If non-GAFF water is present, load desired parameters for that water as well.
    amberparams_section = "\n".join("loadamberparams %s" % (filename) for k, filename in enumerate(frcmod_filenames))

    if solvation == True:
        tleap_commands = TLEAP_TEMPLATE_SOLV % dict(mol2_section=mol2_section, amberparams_section=amberparams_section, buffer_distance=buffer_distance,
                                            box_filename=box_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
    elif solvation == False:
        tleap_commands = TLEAP_TEMPLATE % dict(mol2_section=mol2_section, amberparams_section=amberparams_section,
                                            box_filename=box_filename, prmtop_filename=prmtop_filename, inpcrd_filename=inpcrd_filename)
    print(tleap_commands)

    file_handle = open('tleap_commands', 'w')
    file_handle.writelines(tleap_commands)
    file_handle.close()

    logger.debug('Running tleap in temporary directory.')
    cmd = "tleap -f %s " % file_handle.name
    logger.debug(cmd)

    output = getoutput(cmd)
    logger.debug(output)

    return tleap_commands


# MAIN CODE	
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-i","--infile", type=open, help="Input file that contains names and number of molecules", required=True)
#parser.add_argument('--antibody', type = open)
parser.add_argument("-s", "--solvate", action="store_true", help="If activated, creates a model solvated in explicit water (TIP3P)")
parser.add_argument("-b", "--buffer", type=int, default=10, help="Buffering distance between waterbox and molecule box (default=10 Angstroms)")
parser.add_argument("-#", "--seed", type=int, default=-1, help="If specified, creates a configuration with given seed number (default=random)")
parser.add_argument('-n', '--name', type=str, help = 'Name of input .pdb file containing nanoparticle and antibody.', required = True)
parser.add_argument('-o', '--output', type=str, help = 'Name for output .prmtop and .inpcrd files.', required = True)
args = parser.parse_args()

pair_df = pd.read_csv(args.infile, sep='\t',header=None)
splitpoint = len(pair_df.columns)/2
mol_df = pair_df.loc[:,0:splitpoint-1]
num_df = pair_df.loc[:,splitpoint:]


if args.solvate:
    solvation = True
    buffer_distance = args.buffer
else:
    solvation = False
    buffer_distance = None

seed=args.seed

for i,row in pair_df.iterrows():
    pdb_list = []
    num_list = []
    mol2_filenames = []
    frcmod_filenames = []
    mix_name = ''
    for mol in mol_df.loc[i]:
        mol2_to_pdb(mol+'.mol2')
        pdb_list += [mol+'.pdb']
        mol2_filenames += [mol+'.mol2']
        frcmod_filenames += [mol+'.frcmod']
        mix_name = mix_name + mol + '_'
    for num in num_df.loc[i]:
        num_list += [int(num)]
    
    
    final_name = args.output
    box_filename = args.name


tleap_cmd = build_mixture_prmtop(mol2_filenames, frcmod_filenames, box_filename, final_name+'.prmtop', final_name+'.inpcrd', solvation=solvation, buffer_distance=buffer_distance)

file_handle = open(final_name+".leap.in", 'w')
file_handle.writelines(tleap_cmd)
file_handle.close()


