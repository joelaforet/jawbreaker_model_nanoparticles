"""
@author: Joe Laforet
Created on Mon May 22 12:57:08 2023


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
import hashlib

def map_string_to_color(string):
    # Use the MD5 hashing algorithm to generate a hash value
    hash_object = hashlib.md5(string.encode())
    hash_value = hash_object.hexdigest()

    # Take the first 6 characters of the hash value to represent the color code
    color_code = "#" + hash_value[:6]

    return color_code

def hex_to_rgb(color_code):
    # Remove the '#' symbol if present
    color_code = color_code.lstrip('#')
    
    # Convert hexadecimal to RGB values
    r = int(color_code[0:2], 16)
    g = int(color_code[2:4], 16)
    b = int(color_code[4:6], 16)
    
    return r, g, b

def rgb_to_hex(rgb):
    # Convert RGB values to hexadecimal
    color_code = '#%02x%02x%02x' % rgb
    
    return color_code

def bias_color_towards_white(color_code, bias_factor):
    # Convert color code to RGB values
    r, g, b = hex_to_rgb(color_code)

    # Calculate bias towards white
    r_bias = (255 - r) * bias_factor
    g_bias = (255 - g) * bias_factor
    b_bias = (255 - b) * bias_factor

    # Apply the bias
    r = int(r + r_bias)
    g = int(g + g_bias)
    b = int(b + b_bias)

    # Ensure values are within the valid range
    r = min(max(r, 0), 255)
    g = min(max(g, 0), 255)
    b = min(max(b, 0), 255)

    # Convert RGB back to color code
    biased_color_code = rgb_to_hex((r, g, b))

    return biased_color_code

def bias_color_towards_darker(color_code, bias_factor):
    # Convert color code to RGB values
    r, g, b = hex_to_rgb(color_code)

    # Adjust RGB values towards darker hues
    r = int(r * (1 - bias_factor))
    g = int(g * (1 - bias_factor))
    b = int(b * (1 - bias_factor))

    # Convert RGB back to color code
    biased_color_code = rgb_to_hex((r, g, b))

    return biased_color_code

def bias_dark_color_towards_appealing(color_code, bias_factor):
    # Convert color code to RGB values
    r, g, b = hex_to_rgb(color_code)

    # Calculate the brightness of the color
    brightness = (r + g + b) / 3

    # Calculate the desired target brightness based on the bias factor
    target_brightness = 255 * bias_factor

    # Calculate the adjustment factor for brightness
    adjustment_factor = target_brightness / brightness

    # Adjust the RGB values to steer towards the appealing color
    r = int(r * adjustment_factor)
    g = int(g * adjustment_factor)
    b = int(b * adjustment_factor)

    # Ensure values are within the valid range
    r = min(max(r, 0), 255)
    g = min(max(g, 0), 255)
    b = min(max(b, 0), 255)

    # Convert RGB back to color code
    biased_color_code = rgb_to_hex((r, g, b))

    return biased_color_code



#%% Initialize Argument Parser
parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i","--infile", type=open, help="Input .pdb file that contains simulation trajectory", required=True)
parser.add_argument("-n", "--name", help = "Name of folder to render images to.", required = True)
parser.add_argument("--dr", help = "Residue ID for drug molecule.", required = True)
parser.add_argument("--er", help = "Residue ID for excipient molecule.", required = True)
parser.add_argument("-d", "--delete", action="store_true", default = False, help="If activated, will delete images after rendering video.")
parser.add_argument('-f', '--framerate', type=int, help='Framerate for output movie.', required = True)
args = parser.parse_args()
#print(args.infile)

folder = args.name

#%% Analysis
print("before infile assingnde")
f = args.infile
print("after")
#folder = f.name

def process_pdb_file(input_file, output_file):
    with open(input_file, 'r') as input_pdb:
        lines = input_pdb.readlines()
        
    filtered_lines = [line for line in lines if not line.startswith("ATOM") and not line[17:20] == "HOH"]
    
    with open(output_file, 'w') as output_pdb:
        output_pdb.writelines(filtered_lines)
        
    return output_file

name = f.name.split("_")[2] + '_' + f.name.split("_")[7] + '_' + f.name.split("_")[8]+"_" +f.name.split("_")[-1].split(".")[0]

if not os.path.exists(f"reduced_{name}.pdb"):
    small_file = process_pdb_file(f.name, f"reduced_{name}.pdb")
else:
    small_file = f"reduced_{name}.pdb"

#%% Make and Change to the new folder where we're going to render the images.
curDir = os.getcwd()

if not os.path.exists(folder):
    os.makedirs(folder)
os.chdir(folder)

#%% New folder is made and we are currently inside of it.
## Now we pass the script to pyMol and render the trajectory

def cleanSim(d_resid, d_color, e_resid, e_color):
    """
    Cleans the simulation, colors drug d_color, 
    excipient e_color, and sets Zoom
    
    @param d_resid: input 3 char residue ID of drug
    @param d_color: input 6 char color code hash for drug
    @param e_resid: input 3 char residue ID of excipient
    @param e_color: 6 char color code hash for excipient

    
    Returns
    -------
    None.

    """
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "on")
    cmd.set("orthoscopic", "1")
    cmd.set("depth_cue", "0")
    cmd.set("ray_shadow", "0")
    cmd.remove("solvent")
    cmd.set("specular", 0)
    cmd.show(representation = "spheres", selection = "all")
    cmd.color(f"0x{d_color}",  f"resname {d_resid}")
    cmd.color(f"0x{e_color}",  f"resname {e_resid}")
    cmd.color("yellow", "resname STK")
    cmd.set("sphere_transparency", 0.6, f"not resname {d_resid} and not resname {e_resid}")
    cmd.center(f"resname {d_resid}")
    cmd.orient()
    cmd.zoom("center", "65")
    cmd.smooth("all", 30, 3)


cmd.extend("cleanSim", cleanSim)

####

os.chdir(curDir)

#%% Clean the simulation frames
d_resid = args.dr
e_resid = args.er

d_name = f.name.split("_")[7]
e_name = f.name.split("_")[8]

light_bias_factor = 0.7 #for how light it is
dark_bias_factor = 0.3 #for how dark it is
bright_bias_factor = 0.5 #for how much to brighten the dark so it looks pretty

d_color_code = map_string_to_color(d_name)
d_biased_color_code = bias_color_towards_white(d_color_code, light_bias_factor)[1:]

e_color_code = map_string_to_color(e_name)
dark_biased_color_code = bias_color_towards_darker(e_color_code, dark_bias_factor)
e_biased_color_code = bias_dark_color_towards_appealing(dark_biased_color_code, bright_bias_factor)[1:]

os.system("pymol -cq clean_jawBreaker.py")
cmd.load(small_file)

cleanSim(d_resid, d_biased_color_code, e_resid, e_biased_color_code)

#%% Render the frames to the current working directory
cmd.set("hash_max", "350")

os.chdir(folder)
cmd.mpng(name+ "_", width= 1014, height = 720, mode = 2)
pymol.finish_launching()

framerate = args.framerate
delete = args.delete

os.system("ffmpeg -r {} -i {}_%04d.png -c:v libx264 -pix_fmt yuv420p -y {}.mp4".format(framerate, name, name))
if(delete):
	os.system(f'find . -name "*{name}_*.png" -type f -delete')

os.chdir(curDir)
