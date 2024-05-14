"""
@author: Joe Laforet
Created on Tue Jun 15 12:56:18 2021


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
parser.add_argument('-n', '--outname', help = 'Input name for output plots. Form = Drug_#D_#E', required = True)
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

#%% Iterate through every folder in input directory

# Get a way to go back to the original folder
os.chdir(args.directory)
curDir = os.getcwd()
name = args.outname
drug = name.split('_')[0]

temp = []

###################################################################################

# Start with TotHBperFrame
os.chdir(drug + '_TotHBperFrame')

for x in sorted(glob.glob('*TotBPF.png')):
    temp.append(x)
#print(temp)
#print(len(temp))
img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_TotBPF_FULL.png').format(name))

os.chdir(drug + '_TotHBperFrame')

temp = []

for x in sorted(glob.glob('*TotBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_TotBPF_FULL.png').format(name))

temp = []


##########################################################################

os.chdir(drug + '_TotHBperTime')

for x in sorted(glob.glob('*TotTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_TotTHB_FULL.png').format(name))

os.chdir(drug + '_TotHBperTime')
temp = []

for x in sorted(glob.glob('*TotTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_TotTHB_FULL.png').format(name))

temp = []


############################################################################

os.chdir(drug +'_ddHBperFrame')

for x in sorted(glob.glob('*ddBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_ddBPF_FULL.png').format(name))

os.chdir(drug +'_ddHBperFrame')
temp = []

for x in sorted(glob.glob('*ddBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_ddBPF_FULL.png').format(name))

temp = []


##############################################################################

os.chdir(drug +'_ddHBperTime')

for x in sorted(glob.glob('*ddTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_ddTHB_FULL.png').format(name))

os.chdir(drug +'_ddHBperTime')
temp = []

for x in sorted(glob.glob('*ddTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_ddTHB_FULL.png').format(name))

temp = []


################################################################################

os.chdir(drug +'_deHBperFrame')

for x in sorted(glob.glob('*deBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_deBPF_FULL.png').format(name))

os.chdir(drug +'_deHBperFrame')
temp = []

for x in sorted(glob.glob('*deBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_deBPF_FULL.png').format(name))

temp = []


################################################################################

os.chdir(drug +'_deHBperTime')

for x in sorted(glob.glob('*deTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_deTHB_FULL.png').format(name))

os.chdir(drug +'_deHBperTime')
temp = []

for x in sorted(glob.glob('*deTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_deTHB_FULL.png').format(name))

temp = []


################################################################################

os.chdir(drug +'_eeHBperFrame')

for x in sorted(glob.glob('*eeBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_eeBPF_FULL.png').format(name))

os.chdir(drug +'_eeHBperFrame')
temp = []

for x in sorted(glob.glob('*eeBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_eeBPF_FULL.png').format(name))

temp = []


################################################################################

os.chdir(drug +'_eeHBperTime')

for x in sorted(glob.glob('*eeTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_eeTHB_FULL.png').format(name))

os.chdir(drug +'_eeHBperTime')
temp = []

for x in sorted(glob.glob('*eeTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_eeTHB_FULL.png').format(name))

temp = []


################################################################################

os.chdir(drug +'_edHBperFrame')

for x in sorted(glob.glob('*edBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_edBPF_FULL.png').format(name))

os.chdir(drug +'_edHBperFrame')
temp = []

for x in sorted(glob.glob('*edBPF.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_edBPF_FULL.png').format(name))

temp = []



################################################################################

os.chdir(drug +'_edHBperTime')

for x in sorted(glob.glob('*edTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_edTHB_FULL.png').format(name))

os.chdir(drug +'_edHBperTime')
temp = []

for x in sorted(glob.glob('*edTHB.png')):
    temp.append(x)

img = image_grid(temp, int(args.rows), int(args.columns))
os.chdir(curDir)
img.save(('{}_edTHB_FULL.png').format(name))

temp = []



################################################################################





