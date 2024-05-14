"""
@author: Joe Laforet
Created on Tue May 25 14:51:50 2021


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

#%% Initialize Argument Parser
parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i","--infile", type=open, help="Input .pdb file that contains simulation trajectory", required=True)
parser.add_argument("-n", "--name", help = "Name of outfile.", required = True)
args = parser.parse_args()
#%% Analysis

f = args.infile
nameFile= f.name

t = md.load(nameFile) # Load trajectory file
topology = md.load(nameFile).topology
filetag = nameFile[7:-4]
    
n_frames = len(t) # Get number of frames in simulation
    
table, bonds = topology.to_dataframe()
reses = table['resName'].unique()[0:2]
    
drugres = reses[0]
excipres = reses[1]
    
    
##% Total H-bonds Counters
total_bonds = 0 #Initialize total bond count
frame_totbonds = [] # Initialize total bonds per frame
running_totbonds = [] # Initialize total bonds over time
    
##% Drug-Drug H-bonds Counters
ddbonds = 0 #Initialize total drug-drug bond count
frame_ddbonds = [] # Initialize drug-drug bonds per frame
running_ddbonds = [] #Initialize drug-drug bonds over time
    
##% Drug-Excip H-bonds Counters
debonds = 0 #Initialize total drug-Excip bond count
frame_debonds = [] # Initialize drug-Excip bonds per frame
running_debonds = [] #Initialize drug-Excip bonds over time
    
##% Excip-Drug H-bonds Counters
edbonds = 0 #Initialize total Excip-Drug bond count
frame_edbonds = [] # Initialize Excip-Drug bonds per frame
running_edbonds = [] #Initialize Excip-Drug bonds over time
    
##% Excip-Excip H-Bonds Counters
eebonds = 0 #Initialize total Excip-Excip bond count
frame_eebonds = [] # Initialize Excip-Excip bonds per frame
running_eebonds = [] #Initialize Excip-Excip bonds over time

NN_Bonds = 0
NO_Bonds = 0
ON_Bonds = 0
OO_Bonds = 0

NF_Bonds = 0
OF_Bonds = 0
FN_Bonds = 0
FO_Bonds = 0
FF_Bonds = 0

for i in range(n_frames):
    working = t[i]  #Pick working frame
    hbonds = md.baker_hubbard(working, exclude_water = True)  #Calculate hydrogen bonds in frame
        
    tot_this_frame_bonds = 0
    dd_this_frame_bonds = 0
    de_this_frame_bonds = 0
    ed_this_frame_bonds = 0
    ee_this_frame_bonds = 0
    
    
    for hbond in hbonds:
            
        # Get the two residues we're looking at in the hbond system
        res1 = str(working.topology.atom(hbond[0])) #Get name of first atom in bond pair
        atom1 = working.topology.atom(hbond[0]).name[0]
        res1 = res1[0:3]
        res2 = str(working.topology.atom(hbond[2])) #Get name of second atom in bond pair
        atom2 = working.topology.atom(hbond[2]).name[0]
        res2 = res2[0:3]
        
        
        #If system is Drug-Drug
        if((res1 == drugres) and (res2 == drugres)):
            ddbonds +=1
            dd_this_frame_bonds +=1
            
        #If system is Drug-Excip
        elif((res1 == drugres) and (res2 == excipres)):
            debonds +=1
            de_this_frame_bonds +=1
            if((atom1 == 'N') and (atom2 == 'N')):
                NN_Bonds +=1
            elif((atom1 == 'N') and (atom2 == 'O')):
                NO_Bonds += 1
            elif((atom1 == 'O') and (atom2 == 'N')):
                ON_Bonds += 1
            elif((atom1 == 'O') and (atom2 == 'O')):
                OO_Bonds += 1
                
                
            elif((atom1 == 'F') and (atom2 == 'O')):
                FO_Bonds += 1
            elif((atom1 == 'F') and (atom2 == 'N')):
                FN_Bonds += 1  
            elif((atom1 == 'O') and (atom2 == 'F')):
                OF_Bonds += 1    
            elif((atom1 == 'N') and (atom2 == 'F')):
                NF_Bonds += 1
            elif((atom1 == 'F') and (atom2 == 'F')):
                FF_Bonds += 1
            
        #If system is Excip-Drug
        elif((res1 == excipres) and (res2 == drugres)):
            edbonds +=1
            ed_this_frame_bonds +=1
            
            if((atom1 == 'N') and (atom2 == 'N')):
                NN_Bonds +=1
            elif((atom1 == 'N') and (atom2 == 'O')):
                NO_Bonds += 1
            elif((atom1 == 'O') and (atom2 == 'N')):
                ON_Bonds += 1
            elif((atom1 == 'O') and (atom2 == 'O')):
                OO_Bonds += 1
            
            elif((atom1 == 'F') and (atom2 == 'O')):
                FO_Bonds += 1
            elif((atom1 == 'F') and (atom2 == 'N')):
                FN_Bonds += 1  
            elif((atom1 == 'O') and (atom2 == 'F')):
                OF_Bonds += 1    
            elif((atom1 == 'N') and (atom2 == 'F')):
                NF_Bonds += 1
            elif((atom1 == 'F') and (atom2 == 'F')):
                FF_Bonds += 1                        
            
        #If system is Excip-Excip    
        elif((res1 == excipres) and (res2 == excipres)):
            eebonds +=1
            ee_this_frame_bonds +=1
        #Add to total bond counters
        tot_this_frame_bonds +=1
        total_bonds +=1
        
    # Drug-Drug Arrays
    running_ddbonds.append(ddbonds)
    frame_ddbonds.append(dd_this_frame_bonds)
        
    #Drug-Excip Arrays
    running_debonds.append(debonds)
    frame_debonds.append(de_this_frame_bonds)
        
    #Excip-Drug Arrays
    running_edbonds.append(edbonds)
    frame_edbonds.append(ed_this_frame_bonds)
        
    #Excip-Excip Arrays
    running_eebonds.append(eebonds)
    frame_eebonds.append(ee_this_frame_bonds)
        
    # Total Bond Arrays
    running_totbonds.append(total_bonds)
    frame_totbonds.append(tot_this_frame_bonds)

    
    
    
#######################################################################################################        

drug = str(args.name).split('_')[0]

path = drug
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)
curDir = os.getcwd()


path = drug + '_TotHBperFrame'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

#%% Output Total bonds/frame as a .csv for further analysis
out1_array = np.array(frame_totbonds)
string = args.name
df1 = pd.DataFrame(out1_array)
df1.to_csv(string+'_TotHBperFrame.csv')
os.chdir(curDir)
#%% Output Total Bonds/Time elapsed to .csv

path = drug + '_TotHBperTime'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out2_array = np.array(running_totbonds)
df2 = pd.DataFrame(out2_array)
df2.to_csv(string+'_TotHBperTime.csv')
os.chdir(curDir)

#%% Output Total Bonds to a .csv

path = drug + '_TotHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df3 = pd.DataFrame(data = [total_bonds], index = [string])
df3.to_csv(string + '_TotHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total Drug-Drug bonds/frame as a .csv for further analysis

path = drug + '_ddHBperFrame'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out4_array = np.array(frame_ddbonds)
string = args.name
df4 = pd.DataFrame(out4_array)
df4.to_csv(string+'_ddHBperFrame.csv')
os.chdir(curDir)

#%% Output Total Drug-Drug Bonds/Time elapsed to .csv

path = drug + '_ddHBperTime'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out5_array = np.array(running_ddbonds)
df5 = pd.DataFrame(out5_array)
df5.to_csv(string+'_ddHBperTime.csv')
os.chdir(curDir)

#%% Output Total Drug-Drug Bonds to a .csv

path = drug + '_ddHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df6 = pd.DataFrame(data = [ddbonds], index = [string])
df6.to_csv(string + '_ddHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total Drug-Excip bonds/frame as a .csv for further analysis

path = drug + '_deHBperFrame'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out7_array = np.array(frame_debonds)
string = args.name
df7 = pd.DataFrame(out7_array)
df7.to_csv(string+'_deHBperFrame.csv')
os.chdir(curDir)

#%% Output Total Drug-Excip Bonds/Time elapsed to .csv

path = drug + '_deHBperTime'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out8_array = np.array(running_debonds)
df8 = pd.DataFrame(out8_array)
df8.to_csv(string+'_deHBperTime.csv')
os.chdir(curDir)

#%% Output Total Drug-Excip Bonds to a .csv

path = drug + '_deHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df9 = pd.DataFrame(data = [debonds], index = [string])
df9.to_csv(string + '_deHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total Excip-Drug bonds/frame as a .csv for further analysis

path = drug + '_edHBperFrame'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out10_array = np.array(frame_edbonds)
string = args.name
df10 = pd.DataFrame(out10_array)
df10.to_csv(string+'_edHBperFrame.csv')
os.chdir(curDir)

#%% Output Total Drug-Excip Bonds/Time elapsed to .csv

path = drug + '_edHBperTime'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out11_array = np.array(running_edbonds)
df11 = pd.DataFrame(out11_array)
df11.to_csv(string+'_edHBperTime.csv')
os.chdir(curDir)

#%% Output Total Drug-Excip Bonds to a .csv

path = drug + '_edHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [edbonds], index = [string])
df12.to_csv(string + '_edHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total Excip-Excip bonds/frame as a .csv for further analysis

path = drug + '_eeHBperFrame'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out13_array = np.array(frame_eebonds)
string = args.name
df13 = pd.DataFrame(out13_array)
df13.to_csv(string+'_eeHBperFrame.csv')
os.chdir(curDir)

#%% Output Total Excip-Excip Bonds/Time elapsed to .csv

path = drug + '_eeHBperTime'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

out14_array = np.array(running_eebonds)
df14 = pd.DataFrame(out14_array)
df14.to_csv(string+'_eeHBperTime.csv')
os.chdir(curDir)

#%% Output Total Drug-Excip Bonds to a .csv

path = drug + '_eeHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df15 = pd.DataFrame(data = [eebonds], index = [string])
df15.to_csv(string + '_eeHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total N-N HBonds as a .csv for further analysis

path = drug + '_NNHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [NN_Bonds], index = [string])
df12.to_csv(string + '_NNHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total N-O HBonds as a .csv for further analysis

path = drug + '_NOHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [NO_Bonds], index = [string])
df12.to_csv(string + '_NOHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total O-N HBonds as a .csv for further analysis

path = drug + '_ONHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [ON_Bonds], index = [string])
df12.to_csv(string + '_ONHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total O-O HBonds as a .csv for further analysis

path = drug + '_OOHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [OO_Bonds], index = [string])
df12.to_csv(string + '_OOHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total F-O HBonds as a .csv for further analysis

path = drug + '_FOHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [FO_Bonds], index = [string])
df12.to_csv(string + '_FOHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total F-N HBonds as a .csv for further analysis

path = drug + '_FNHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [FN_Bonds], index = [string])
df12.to_csv(string + '_FNHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total O-F HBonds as a .csv for further analysis

path = drug + '_OFHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [OF_Bonds], index = [string])
df12.to_csv(string + '_OFHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total N-F HBonds as a .csv for further analysis

path = drug + '_NFHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [NF_Bonds], index = [string])
df12.to_csv(string + '_NFHBonds.csv')
os.chdir(curDir)

#######################################################################################################

#%% Output Total F-F HBonds as a .csv for further analysis

path = drug + '_FFHBonds'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

df12 = pd.DataFrame(data = [FF_Bonds], index = [string])
df12.to_csv(string + '_FFHBonds.csv')
os.chdir(curDir)

#########################################################################################

