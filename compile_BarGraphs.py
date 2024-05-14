"""
@author: Joe Laforet
Created on Wed Jun 16 13:22:18 2021


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""

#%% Import Modules
import argparse
import PIL, os, glob
from PIL import Image
from math import ceil, floor

#%% Initialize Argparse

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument('-d', '--directory', help = 'Input directory where files are stored/ plots will be stored in.', required = True)
parser.add_argument('-n', '--name', help = 'Input name for output plots. Form = Drug_#D_#E', required = True)
parser.add_argument('-r', '--rows', help = 'Input how many rows in compiled figures.', required = True)
parser.add_argument('-c', '--columns', help = 'Input how many columns in compiled figures.', required = True)

args = parser.parse_args()

#%% Image making function

def image_grid(imgs, rows, cols):
    #print(len(imgs))
    #print(rows*cols)
    #assert len(imgs) == (rows*cols)
    
    w, h = Image.open(imgs[0]).size
    grid = Image.new('RGB', size=(cols*w, rows*h))
    grid_w, grid_h = grid.size
    
    for i, img in enumerate(imgs):
        img = Image.open(img)
        grid.paste(img, box=(i%cols*w, i//cols*h))
    return grid

#%%

path = args.directory
name = args.name

os.chdir(path)

curDir = os.getcwd()
os.chdir('BarGraphs')

temp = []

##################################################################

#Start with Total H-Bond Graphs

for x in sorted(glob.glob('*TotalHBonds*.png')):
    temp.append(x)
#print(len(temp))
#print(int(args.rows) * int(args.columns))
img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_TotHBonds_FULL.png').format(name))

temp = []
os.chdir('BarGraphs')
##################################################################

#D-D H-Bond Graphs

for x in sorted(glob.glob('*ddHBonds*.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_ddHBonds_FULL.png').format(name))

temp = []
os.chdir('BarGraphs')
##################################################################

#D-E H-Bond Graphs

for x in sorted(glob.glob('*deHBonds*.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_deHBonds_FULL.png').format(name))

temp = []
os.chdir('BarGraphs')
##################################################################

#E-D H-Bond Graphs

for x in sorted(glob.glob('*edHBonds*.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_edHBonds_FULL.png').format(name))

temp = []
os.chdir('BarGraphs')
##################################################################

#E-E H-Bond Graphs

for x in sorted(glob.glob('*eeHBonds*.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_eeHBonds_FULL.png').format(name))

temp = []
os.chdir('BarGraphs')
##################################################################

#N-N H-Bond Graphs

for x in sorted(glob.glob('*NNHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_NNHBonds_FULL.png').format(name))

#temp = []

##################################################################

#N-O H-Bond Graphs

for x in sorted(glob.glob('*NOHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_NOHBonds_FULL.png').format(name))

#temp = []

##################################################################

#O-N H-Bond Graphs

for x in sorted(glob.glob('*ONHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_ONHBonds_FULL.png').format(name))

#temp = []

##################################################################

#O-O H-Bond Graphs

for x in sorted(glob.glob('*OOHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_OOHBonds_FULL.png').format(name))

#temp = []

##################################################################

#F-O H-Bond Graphs

for x in sorted(glob.glob('*FOHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_FOHBonds_FULL.png').format(name))

#temp = []

##################################################################

#F-N H-Bond Graphs

for x in sorted(glob.glob('*FNHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_FNHBonds_FULL.png').format(name))

#temp = []

##################################################################

#O-F H-Bond Graphs

for x in sorted(glob.glob('*OFHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_OFHBonds_FULL.png').format(name))

#temp = []

##################################################################

#N-F H-Bond Graphs

for x in sorted(glob.glob('*NFHBonds*.png')):
    temp.append(x)

#img = image_grid(temp, int(args.rows), int(args.columns))
#os.chdir(curDir)
#img.save(('{}_NFHBonds_FULL.png').format(name))

#temp = []

##################################################################

#F-F H-Bond Graphs

for x in sorted(glob.glob('*FFHBonds*.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_Atomic_HBonds_FULL.png').format(name))

temp = []

###################################################################
