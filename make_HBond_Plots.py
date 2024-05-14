"""
@author: Joe Laforet
Created on Thu May 27 15:22:25 2021


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""

#%% Import Modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

#%% Initialize Argparse

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument('-n', '--outname', help = "Input name for output files.", required = True)
parser.add_argument('-d', '--directory', help = "Input directory where data is stored.", required = True)

args = parser.parse_args()

########################################################################################

os.chdir(args.directory)
curDir = os.getcwd()
name = args.outname
drug = name.split('_')[0]

#%% Setting up Total HBonds/Frame Dataframe
os.chdir(drug + '_TotHBperFrame')

frames = args.outname + '_TotHBperFrame.csv'

Frame_bonds = pd.read_csv(frames)
Frame_bonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for HBonds/Frame
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total HBonds/Frame for {}'.format(args.outname))
ax1.grid(True)
ax1.plot(Frame_bonds['Frame'], Frame_bonds['Bonds'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_TotBPF.png').format(args.outname))

os.chdir(curDir)
#%% Setting up Total HBonds 
os.chdir(drug + '_TotHBperTime')
total = args.outname + '_TotHBperTime.csv'

totalbonds = pd.read_csv(total)
totalbonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for Total HBonds
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig2.add_subplot(1,1,1)
ax2.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total Hydrogen Bonds for {}'.format(args.outname))
ax2.grid(True)
ax2.plot(totalbonds['Frame'], totalbonds['Bonds'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_TotTHB.png').format(args.outname))
os.chdir(curDir)
########################################################################################

#%% Setting up Drug-Drug HBonds/Frame Dataframe
os.chdir(drug +'_ddHBperFrame')
frames = args.outname + '_ddHBperFrame.csv'

Frame_bonds = pd.read_csv(frames)
Frame_bonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for HBonds/Frame
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Drug-Drug HBonds/Frame for {}'.format(args.outname))
ax1.grid(True)
ax1.plot(Frame_bonds['Frame'], Frame_bonds['Bonds'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_ddBPF.png').format(args.outname))
os.chdir(curDir)

#%% Setting up Total HBonds 
os.chdir(drug +'_ddHBperTime')
total = args.outname + '_ddHBperTime.csv'

totalbonds = pd.read_csv(total)
totalbonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for Total HBonds
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig2.add_subplot(1,1,1)
ax2.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total Drug-Drug Hydrogen Bonds for {}'.format(args.outname))
ax2.grid(True)
ax2.plot(totalbonds['Frame'], totalbonds['Bonds'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_ddTHB.png').format(args.outname))
os.chdir(curDir)
########################################################################################

#%% Setting up Drug-Excipient HBonds/Frame Dataframe
os.chdir(drug +'_deHBperFrame')
frames = args.outname + '_deHBperFrame.csv'

Frame_bonds = pd.read_csv(frames)
Frame_bonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for HBonds/Frame
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Drug-Excipient HBonds/Frame for {}'.format(args.outname))
ax1.grid(True)
ax1.plot(Frame_bonds['Frame'], Frame_bonds['Bonds'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_deBPF.png').format(args.outname))
os.chdir(curDir)

#%% Setting up Total HBonds 
os.chdir(drug +'_deHBperTime')
total = args.outname + '_deHBperTime.csv'

totalbonds = pd.read_csv(total)
totalbonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for Total HBonds
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig2.add_subplot(1,1,1)
ax2.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total Drug-Excipient Hydrogen Bonds for {}'.format(args.outname))
ax2.grid(True)
ax2.plot(totalbonds['Frame'], totalbonds['Bonds'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_deTHB.png').format(args.outname))
os.chdir(curDir)
########################################################################################

#%% Setting up Excipient-Excipient HBonds/Frame Dataframe
os.chdir(drug +'_eeHBperFrame')
frames = args.outname + '_eeHBperFrame.csv'

Frame_bonds = pd.read_csv(frames)
Frame_bonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for HBonds/Frame
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Excipient-Excipient HBonds/Frame for {}'.format(args.outname))
ax1.grid(True)
ax1.plot(Frame_bonds['Frame'], Frame_bonds['Bonds'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_eeBPF.png').format(args.outname))
os.chdir(curDir)

#%% Setting up Total HBonds 
os.chdir(drug +'_eeHBperTime')
total = args.outname + '_eeHBperTime.csv'

totalbonds = pd.read_csv(total)
totalbonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for Total HBonds
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig2.add_subplot(1,1,1)
ax2.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total Excipient-Excipient Hydrogen Bonds for {}'.format(args.outname))
ax2.grid(True)
ax2.plot(totalbonds['Frame'], totalbonds['Bonds'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_eeTHB.png').format(args.outname))
os.chdir(curDir)
########################################################################################

#%% Setting up Excipient-Drug HBonds/Frame Dataframe
os.chdir(drug +'_edHBperFrame')
frames = args.outname + '_edHBperFrame.csv'

Frame_bonds = pd.read_csv(frames)
Frame_bonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for HBonds/Frame
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Excipient-Drug HBonds/Frame for {}'.format(args.outname))
ax1.grid(True)
ax1.plot(Frame_bonds['Frame'], Frame_bonds['Bonds'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_edBPF.png').format(args.outname))
os.chdir(curDir)

#%% Setting up Total HBonds 
os.chdir(drug +'_edHBperTime')
total = args.outname + '_edHBperTime.csv'

totalbonds = pd.read_csv(total)
totalbonds.rename(columns = {'Unnamed: 0': 'Frame', '0': 'Bonds'}, inplace = True)

#%% Make Figure for Total HBonds
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig2.add_subplot(1,1,1)
ax2.set(xlabel = 'Frame Number', ylabel = 'HBonds formed', title = 'Total Excipient-Drug Hydrogen Bonds for {}'.format(args.outname))
ax2.grid(True)
ax2.plot(totalbonds['Frame'], totalbonds['Bonds'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_edTHB.png').format(args.outname))
os.chdir(curDir)
###########################################################################################
















