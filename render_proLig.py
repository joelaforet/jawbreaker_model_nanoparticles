"""
@author: Joe Laforet
Created on Sat June 18 13:36:08 2022


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""



#%% Import Modules
import mdtraj as md
import os
import argparse
import csv
import numpy as np
import pandas as pd
from pymol import cmd,stored

#%% Initialize Argument Parser
parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i","--infile", type=open, help="Input .pdb file that contains simulation trajectory", required=True)
parser.add_argument("-n", "--name", help = "Name of folder to render images to.", required = True)
parser.add_argument("-r", "--resid", help = "Residue ID for ligand.", required = True)
parser.add_argument("-d", "--del", action="store_true", default = False, help="If activated, will delete images after rendering video.")
parser.add_argument('-f', '--framerate', type=int, help='Framerate for output movie.', required = True)
args = parser.parse_args()
#print(args.infile)

folder = args.name

#%% Analysis

f = args.infile

#%% Make and Change to the new folder where we're going to render the images.
curDir = os.getcwd()

if not os.path.exists(folder):
    os.makedirs(folder)
os.chdir(folder)

#%% New folder is made and we are currently inside of it.
## Now we pass the script to pyMol and render the trajectory

def cleanSim(resid):
    """
    Cleans the simulation, colors protein gray, 
    ligand by atomID, and sets Zoom
    
    @param resNames: input residue ID of ligand
    
    Returns
    -------
    None.

    """

    tag = "resname {}".format(resid)
 
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "on")
    cmd.remove("solvent")
    cmd.show(representation = "surface", selection = "not resname {}".format(resid))
    cmd.color("gray90",  "all")
    cmd.show(representation = "licorice", selection = tag)
    cmd.color("atomic" , selection = tag)

    cmd.intra_fit("(name ca)")
    cmd.set("transparency", "0.5", "not resname MCB")
    cmd.center(tag)
    cmd.orient()
    cmd.zoom("center", "30")
    cmd.turn("y", "180") 

cmd.extend("cleanSim", cleanSim)



####

os.chdir(curDir)


os.system("pymol -cq clean.py")
cmd.load(f.name)
#%% Clean the simulation frames
resid = args.resid

cleanSim(resid)

#%% Render the frames to the current working directory
cmd.set("hash_max", "350")
os.chdir(folder)
cmd.mpng(folder+ "_", width= 1014, height = 720, mode = 2)
pymol.finish_launching()

framerate = args.framerate
name = args.name
delete = args.del

os.system("ffmpeg -r {} -i {}_%04d.png -c:v libx264 -pix_fmt yuv420p -y {}.mp4".format(framerate, name, name))
if(delete):
	os.system('find . -name "*.png" -type f -delete')

os.chdir(curDir)
