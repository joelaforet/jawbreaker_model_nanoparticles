"""
pam_pro_lig.py
by Joe Laforet Jr.
jrl78@duke.edu

"""
'''
Generate AMBER-recognized input files for a protein-ligand system.
usage: pam_pro_lig.py [-h] --complex COMPLEX --ligand LIGAND -r RESID [-s] [-b BUFFER]
	optional arguments:
	  -h, --help            show this help message and exit
	  -c, --complex         Input file of the protein-ligand complex
	  -l, --ligand          Input file of the ligand
          -r, --resID           Input residue ID code associated with ligand
	  -s, --solvate         If activated, creates a model solvated in explicit water (TIP3P)
	  -b BUFFER, --buffer BUFFER
				Buffering distance between waterbox and molecule box
				(default=10 Angstroms)
'''

import pandas as pd
import numpy as np
import os
from distutils.spawn import find_executable
from pymol import cmd,stored
import sys
import argparse

#############################################################################

def protonate(molecule):
    """
    Protonates input molecule.

    @param ligand: input name of molecule to be protonated

    Returns
    -------
    Protonated .pdb file of molecule
    """
    
    cmd.load(molecule + ".pdb")
    cmd.h_add()
    cmd.save(molecule + "_h.pdb")

##############################################################################

def make_mol2(ligand):
    """
    Generates .mol2 and .frcmod files for ligand.

    @param ligand: input protonated ligand .pdb file

    Returns
    -------
    .mol2 and .frcmod files for input ligand.
    """
    print("Running Antechamber...")
    print("###############################################################################################################")
    os.system("antechamber -i {}_h.pdb -fi pdb -o {}.mol2 -fo mol2 -c bcc -s 2".format(ligand, ligand))
    if( not os.path.exists("{}.mol2".format(ligand))):
        print("Antechamber failed, odd number of electrons detected.")
        print("Running ObabeWYF parameter for cation-pi interactionsl...")
        os.system("obabel -ipdb {}_h.pdb -omol2 -O {}.mol2".format(ligand,ligand))
        print("Running Antechamber to convert atom types...")
        print("###############################################################################################################")
        os.system("antechamber -i {}.mol2 -fi mol2 -o {}.mol2 -fo mol2".format(ligand, ligand))
    print("###############################################################################################################")
    print("Generating .frcmod file...")
    print("Running parmchk2...")
    os.system("parmchk2 -i {}.mol2 -f mol2 -o {}.frcmod -s gaff2 -a Y".format(ligand, ligand))
    print("Parmchk2 run successful!")
    print("Ligand successfully parametrized!")

###############################################################################

def pam_ligand(resID, ligand):

    """
    Generates .lib file for ligand.

    @param ligand: input name of ligand

    Returns
    -------
    .lib file for input ligand.
    """

    template = """
    source leaprc.gaff2
    {} = loadmol2 {}.mol2
    loadamberparams {}.frcmod
    saveoff {} {}.lib
    quit
    """.format(resID, ligand, ligand, resID, ligand)

    file_handle = open(ligand + ".leap.in", 'w')
    file_handle.writelines(template)
    file_handle.close()

    cmd = "tleap -f %s " % file_handle.name

    os.system(cmd)

#############################################################################

def pam_complex(resID, ligand, combined, solvate = True, buffer = 10.0):
    """
    Generates .prmtop and .inpcrd files for protein-ligand complex.

    @param resID: input residue code for ligand
    @param ligand: input name of ligand
    @param solvate: (boolean) solvate in explicit water
    @param buffer: specify degree of padding for water box (angstroms)

    Returns
    -------
    .prmtop and .inpcrd files for protein-ligand system.
    """
    if solvate:
        template = """
        source leaprc.gaff2
        source leaprc.protein.ff19SB
        source leaprc.water.tip3p
        loadamberparams {}.frcmod
        loadoff {}.lib
        box = loadpdb {}.pdb
        addions box Na+ 0
        addions box Cl- 0
        solvateoct box TIP3PBOX {:d}
        saveamberparm box {}_solv.prmtop {}_solv.inpcrd
        quit
        """.format(ligand, ligand, combined, buffer, combined, combined)
    else:
        template = """
        source leaprc.gaff2
        source leaprc.protein.ff19SB
        source leaprc.water.tip3p
        loadamberparams {}.frcmod
        loadoff {}.lib
        box = loadpdb {}.pdb
        addions box Na+ 0
        addions box Cl- 0
        saveamberparm box {}_solv.prmtop {}_solv.inpcrd
        quit
        """.format(ligand, ligand, combined, combined, combined)
    file_handle = open(combined + ".leap.in", 'w')
    file_handle.writelines(template)
    file_handle.close()
    cmd = "tleap -f %s " % file_handle.name
    print("Running Leap...")
    print("###############################################################################################################")
    os.system(cmd)    
    print("Leap Run Successful!")

#################################################################################################################################################
parser = argparse.ArgumentParser(add_help=True)
parser.add_argument("-c", "--complex", type=str, help="Input name of the protein-ligand complex", required=True)
parser.add_argument("-l", "--ligand", type=str, help="Input name of the ligand", required=True)
parser.add_argument("-r","--resID", type=str, help="Input residue ID code associated with ligand", required=True)
parser.add_argument("-s", "--solvate", type=bool, default = True, help="If activated, creates a model solvated in explicit water (TIP3P)")
parser.add_argument("-b", "--buffer", type=int, default=10, help="Buffering distance between waterbox and molecule box (default=10 Angstroms)")
args = parser.parse_args()
################################################################################################################################################

comp = args.complex
lig = args.ligand
res = args.resID
solv = args.solvate
buff = args.buffer

print("Beginning Protonation...")

protonate(lig)

print("Protonation Successful!")

print("Beginning Ligand Parametrization...")

make_mol2(lig)

print("###############################################################################################################")
print("Running Leap for the Ligand...")
print("###############################################################################################################")
pam_ligand(res, lig)

print("Leap Run Successful!")

print("###############################################################################################################")
print("Protonating Complex...")

cmd.load(comp + ".pdb")
cmd.h_add("resname {}".format(res))
cmd.save(comp + "_h.pdb")

comp = comp + "_h"
print("Protonation Successful!")
print("###############################################################################################################")
print("Running Leap for the Complex...")
print("###############################################################################################################")
pam_complex(res, lig, comp, solv, buff)

print("Leap Run Successful!")
print("###############################################################################################################")
print("Job's Done!")
