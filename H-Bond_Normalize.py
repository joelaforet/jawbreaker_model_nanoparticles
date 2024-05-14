"""
@author: Joe Laforet
Created on Wed Jun  2 13:30:34 2021


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""

#%% This script generates a .tsv of unique tags for input drug/excipient pair combinations. 
# The .tsv can then be appended to the larger dataframe of H-Bond analysis.

#%% Import Modules
import pandas as pd
import numpy as np
import rdkit.Chem.Lipinski
from rdkit import Chem
import argparse
import glob
import os

#%% Initialize Argparse
parser = argparse.ArgumentParser(add_help = True)
parser.add_argument('--inmols', help = 'Input .tsv of Drug/Excipient and their SMILES.', required = True)
parser.add_argument('--incomp', help = 'Input .tsv of Drug/Excipient pair and pair composition.', required = True)
parser.add_argument('-n', '--name', help = 'Input Name for the output dataframe.', required = True)
parser.add_argument('-d', '--directory', help = "Input directory where data is stored.", required = True)
args = parser.parse_args()

#%% 
mols = args.inmols
comp = args.incomp

# Make and clean the mols DataFrame
molsdf = pd.read_csv(mols, delimiter = '\t', header = None)
molsdf.rename(columns = {0:'Molecule', 1: 'SMILES'}, inplace = True)

# Make and clean the pairs DataFrame
pairsdf = pd.read_csv(comp, delimiter = '\t', header = None)
pairsdf.rename(columns = {0: 'Drug', 1: 'Excipient', 2: 'Num Drug', 3: 'Num Excipient'}, inplace = True)
pairsdf.sort_values(['Drug', 'Excipient'], inplace = True)

# Make the tag Column in the pairs DataFrame
pairsdf['Num Drug'] = pairsdf['Num Drug'].apply(str)
pairsdf['Num Excipient'] = pairsdf['Num Excipient'].apply(str)
pairsdf['Tag'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'] + '_' +pairsdf['Num Excipient']

# Change the Num Excip/ Num Drug Columns back to integers so we can do math
pairsdf['Num Drug'] = pairsdf['Num Drug'].apply(int)
pairsdf['Num Excipient'] = pairsdf['Num Excipient'].apply(int)


# Make the Drug SMILES column from the mols DataFrame
running = []
for x in set(pairsdf.Drug):
    if x in set(molsdf['Molecule']):
        smiles = (molsdf.loc[molsdf['Molecule'] == x]['SMILES']).to_list()[0]
        for y in range(sum(pairsdf['Drug'] == x)):
            running.append(smiles)
pairsdf['Drug SMILES'] = running

# Make the excipient SMILES column from the mols DataFrame
running = []
for x in set(pairsdf.Drug):
    for y in (pairsdf[pairsdf['Drug']==x]).Excipient:
        smiles = (molsdf.loc[molsdf['Molecule'] == y]['SMILES']).to_list()[0]
        running.append(smiles)
pairsdf['Excipient SMILES'] = running

# Fix the Index column
pairsdf.reset_index(inplace = True)
pairsdf.drop(columns = 'index', inplace = True)

# Make the H-Bond Acceptor/Donor Site columns
pairsdf['Drug H-Acceptors'] = pairsdf.apply(lambda row: Chem.Lipinski.NumHAcceptors(Chem.MolFromSmiles(row['Drug SMILES'])), axis = 1)
pairsdf['Excipient H-Acceptors'] = pairsdf.apply(lambda row: Chem.Lipinski.NumHAcceptors(Chem.MolFromSmiles(row['Excipient SMILES'])), axis = 1)
pairsdf['Drug H-Donors'] = pairsdf.apply(lambda row: Chem.Lipinski.NumHDonors(Chem.MolFromSmiles(row['Drug SMILES'])), axis = 1)
pairsdf['Excipient H-Donors'] = pairsdf.apply(lambda row: Chem.Lipinski.NumHDonors(Chem.MolFromSmiles(row['Drug SMILES'])), axis = 1)

# Make the different combinations of numerical H-Bond Acceptor/Donor site information
pairsdf['Sum Acceptors'] = pairsdf['Drug H-Acceptors'] + pairsdf['Excipient H-Acceptors']
pairsdf['Sum Donors'] = pairsdf['Drug H-Donors'] + pairsdf['Excipient H-Donors']

pairsdf['Sum Excip A/D'] = pairsdf['Excipient H-Acceptors'] + pairsdf['Excipient H-Donors']
pairsdf['Sum Drug A/D'] = pairsdf['Drug H-Acceptors'] + pairsdf['Drug H-Donors']



# Get the total atom count based on SMILES strings for drug/excip system
def atom_count(row):
    drug = Chem.MolFromSmiles(row['Drug SMILES'])
    excip = Chem.MolFromSmiles(row['Excipient SMILES'])
    
    return(drug.GetNumAtoms() + excip.GetNumAtoms())
    
pairsdf['Total Atom Num'] = pairsdf.apply(atom_count, axis = 1)

#################################################################################################################

os.chdir(args.directory)
curDir = os.getcwd()
name = args.name
drug = name.split('_')[0]

path = drug + '_TotHBonds'
os.chdir(path)

#%% Add the total H-Bond Count

TotalBondDF = [x for x in glob.glob('*TotHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df2 = pd.read_csv(x)
    dflist.append(df2)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
pairsdf['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)
tempdf.sort_values(by = 'Temp', inplace = True)

bigdf1 = pd.merge(pairsdf, tempdf, on = 'Temp', how = 'left')
bigdf1.drop(columns = ['Temp'], inplace = True)
bigdf1.rename(columns = {'0': 'Total H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total Drug-Drug H-Bond Count

path = drug + '_ddHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*ddHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf1['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf2 = pd.merge(bigdf1, tempdf, on = 'Temp', how = 'left')
bigdf2.drop(columns = ['Temp'], inplace = True)
bigdf2.rename(columns = {'0': 'Total Drug-Drug H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total Drug-Excipient H-Bond Count

path = drug + '_deHBonds'
os.chdir(path)

TotalBondDF = [x for x in glob.glob('*deHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df1 = pd.read_csv(x)
    dflist.append(df1)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf2['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf3 = pd.merge(bigdf2, tempdf, on = 'Temp', how = 'left')
bigdf3.drop(columns = ['Temp'], inplace = True)
bigdf3.rename(columns = {'0': 'Total Drug-Excipient H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total Excipient-Drug H-Bond Count

path = drug + '_edHBonds'
os.chdir(path)

TotalBondDF = [x for x in glob.glob('*edHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df1 = pd.read_csv(x)
    dflist.append(df1)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf3['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf4 = pd.merge(bigdf3, tempdf, on = 'Temp', how = 'left')
bigdf4.drop(columns = ['Temp'], inplace = True)
bigdf4.rename(columns = {'0': 'Total Excipient-Drug H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total Excipient-Excipient H-Bond Count

path = drug + '_eeHBonds'
os.chdir(path)

TotalBondDF = [x for x in glob.glob('*eeHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df1 = pd.read_csv(x)
    dflist.append(df1)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf4['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf5 = pd.merge(bigdf4, tempdf, on = 'Temp', how = 'left')
bigdf5.drop(columns = ['Temp'], inplace = True)
bigdf5.rename(columns = {'0': 'Total Excipient-Excipient H-Bonds'}, inplace = True)

os.chdir(curDir)
#bigdf5.drop(columns = ['Unnamed: 0'], inplace = True)
#################################################################################################################

#%% Add the total N-N H-Bond Count

path = drug + '_NNHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*NNHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df1 = pd.read_csv(x)
    dflist.append(df1)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf5['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf6 = pd.merge(bigdf5, tempdf, on = 'Temp', how = 'left')
bigdf6.drop(columns = ['Temp'], inplace = True)
#print(bigdf6.columns)
bigdf6.rename(columns = {'0': 'Total N-N H-Bonds'}, inplace = True)
#print(bigdf6.columns)

os.chdir(curDir)

#################################################################################################################

#%% Add the total N-O H-Bond Count

path = drug + '_NOHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*NOHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df1 = pd.read_csv(x)
    dflist.append(df1)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf6['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf7 = pd.merge(bigdf6, tempdf, on = 'Temp', how = 'left')
bigdf7.drop(columns = ['Temp'], inplace = True)
bigdf7.rename(columns = {'0': 'Total N-O H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total O-N H-Bond Count

path = drug + '_ONHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*ONHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf7['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf8 = pd.merge(bigdf7, tempdf, on = 'Temp', how = 'left')
bigdf8.drop(columns = ['Temp'], inplace = True)
bigdf8.rename(columns = {'0': 'Total O-N H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total O-O H-Bond Count

path = drug + '_OOHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*OOHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf8['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf9 = pd.merge(bigdf8, tempdf, on = 'Temp', how = 'left')
bigdf9.drop(columns = ['Temp'], inplace = True)
bigdf9.rename(columns = {'0': 'Total O-O H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total F-O H-Bond Count

path = drug + '_FOHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*FOHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf9['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf10 = pd.merge(bigdf9, tempdf, on = 'Temp', how = 'left')
bigdf10.drop(columns = ['Temp'], inplace = True)
bigdf10.rename(columns = {'0': 'Total F-O H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total F-N H-Bond Count

path = drug + '_FNHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*FNHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf10['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf11 = pd.merge(bigdf10, tempdf, on = 'Temp', how = 'left')
bigdf11.drop(columns = ['Temp'], inplace = True)
bigdf11.rename(columns = {'0': 'Total F-N H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################

#%% Add the total O-F H-Bond Count

path = drug + '_OFHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*OFHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf11['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf12 = pd.merge(bigdf11, tempdf, on = 'Temp', how = 'left')
bigdf12.drop(columns = ['Temp'], inplace = True)
bigdf12.rename(columns = {'0': 'Total O-F H-Bonds'}, inplace = True)
#print('Test 1:')
#print(bigdf12['Total O-F H-Bonds'])
#print(bigdf12.columns)

os.chdir(curDir)

#################################################################################################################

#%% Add the total N-F H-Bond Count

path = drug + '_NFHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*NFHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf12['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf13 = pd.merge(bigdf12, tempdf, on = 'Temp', how = 'left')
bigdf13.drop(columns = ['Temp'], inplace = True)
bigdf13.rename(columns = {'0': 'Total N-F H-Bonds'}, inplace = True)
#print(bigdf13.columns)

os.chdir(curDir)

#################################################################################################################

#%% Add the total F-F H-Bond Count

path = drug + '_FFHBonds'
os.chdir(path)


TotalBondDF = [x for x in glob.glob('*FFHBonds.csv')]
dflist = []

for x in TotalBondDF:
    df3 = pd.read_csv(x)
    dflist.append(df3)
    
tempdf = pd.concat(dflist)

tempdf.rename(columns = {'Unnamed: 0' :'Temp'}, inplace = True)
bigdf13['Temp'] = pairsdf['Drug'] + '_' + pairsdf['Excipient'] + '_' + pairsdf['Num Drug'].apply(str) + '_' + pairsdf['Num Excipient'].apply(str)

bigdf14 = pd.merge(bigdf13, tempdf, on = 'Temp', how = 'left')
bigdf14.drop(columns = ['Temp'], inplace = True)
bigdf14.rename(columns = {'0': 'Total F-F H-Bonds'}, inplace = True)

os.chdir(curDir)

#################################################################################################################


bigdf = bigdf14
#print(bigdf.columns)
#print(bigdf['Total O-F H-Bonds'])

#%% Make the Normalized Total H-Bond Columns

#Excipient Hydrogen Acceptors, Excipient Hydrogen Donors, Total Excipient Acceptors/Donors
bigdf['Norm Tot H-Bonds: ExHA'] = bigdf['Total H-Bonds']/ bigdf['Excipient H-Acceptors']
bigdf['Norm Tot H-Bonds: ExHD'] = bigdf['Total H-Bonds']/ bigdf['Excipient H-Donors']
bigdf['Norm Tot H-Bonds: TExAD'] = bigdf['Total H-Bonds']/ bigdf['Sum Excip A/D']

#Drug Hydrogen Acceptors, Drug Hydrogen Donors, Total Drug Acceptors/Donors
bigdf['Norm Tot H-Bonds: DrHA'] = bigdf['Total H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm Tot H-Bonds: DrHD'] = bigdf['Total H-Bonds']/ bigdf['Drug H-Donors']
bigdf['Norm Tot H-Bonds: TDrAD'] = bigdf['Total H-Bonds']/ bigdf['Sum Drug A/D']

#Total Hydrogen Acceptors, Total Hydrogen Donors
bigdf['Norm Tot H-Bonds: THA'] = bigdf['Total H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm Tot H-Bonds: THD'] = bigdf['Total H-Bonds']/ bigdf['Drug H-Donors']

###########################################################################################

#%% Make the Normalized Total Drug-Drug H-Bond Columns

#Excipient Hydrogen Acceptors, Excipient Hydrogen Donors, Total Excipient Acceptors/Donors
bigdf['Norm D-D H-Bonds: ExHA'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Excipient H-Acceptors']
bigdf['Norm D-D H-Bonds: ExHD'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Excipient H-Donors']
bigdf['Norm D-D H-Bonds: TExAD'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Sum Excip A/D']

#Drug Hydrogen Acceptors, Drug Hydrogen Donors, Total Drug Acceptors/Donors
bigdf['Norm D-D H-Bonds: DrHA'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm D-D H-Bonds: DrHD'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Drug H-Donors']
bigdf['Norm D-D H-Bonds: TDrAD'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Sum Drug A/D']

#Total Hydrogen Acceptors, Total Hydrogen Donors
bigdf['Norm D-D H-Bonds: THA'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm D-D H-Bonds: THD'] = bigdf['Total Drug-Drug H-Bonds']/ bigdf['Drug H-Donors']

###########################################################################################

#%% Make the Normalized Total Drug-Excipient H-Bond Columns

bigdf['Norm D-E H-Bonds: ExHA'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Excipient H-Acceptors']
bigdf['Norm D-E H-Bonds: ExHD'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Excipient H-Donors']
bigdf['Norm D-E H-Bonds: TExAD'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Sum Excip A/D']

#Drug Hydrogen Acceptors, Drug Hydrogen Donors, Total Drug Acceptors/Donors
bigdf['Norm D-E H-Bonds: DrHA'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm D-E H-Bonds: DrHD'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Drug H-Donors']
bigdf['Norm D-E H-Bonds: TDrAD'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Sum Drug A/D']

#Total Hydrogen Acceptors, Total Hydrogen Donors
bigdf['Norm D-E H-Bonds: THA'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm D-E H-Bonds: THD'] = bigdf['Total Drug-Excipient H-Bonds']/ bigdf['Drug H-Donors']

###########################################################################################

#%% Make the Normalized Total Excipient-Drug H-Bond Columns

bigdf['Norm E-D H-Bonds: ExHA'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Excipient H-Acceptors']
bigdf['Norm E-D H-Bonds: ExHD'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Excipient H-Donors']
bigdf['Norm E-D H-Bonds: TExAD'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Sum Excip A/D']

#Drug Hydrogen Acceptors, Drug Hydrogen Donors, Total Drug Acceptors/Donors
bigdf['Norm E-D H-Bonds: DrHA'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm E-D H-Bonds: DrHD'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Drug H-Donors']
bigdf['Norm E-D H-Bonds: TDrAD'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Sum Drug A/D']

#Total Hydrogen Acceptors, Total Hydrogen Donors
bigdf['Norm E-D H-Bonds: THA'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm E-D H-Bonds: THD'] = bigdf['Total Excipient-Drug H-Bonds']/ bigdf['Drug H-Donors']

###########################################################################################

#%% Make the Normalized Total Excipient-Excipient H-Bond Columns

bigdf['Norm E-E H-Bonds: ExHA'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Excipient H-Acceptors']
bigdf['Norm E-E H-Bonds: ExHD'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Excipient H-Donors']
bigdf['Norm E-E H-Bonds: TExAD'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Sum Excip A/D']

#Drug Hydrogen Acceptors, Drug Hydrogen Donors, Total Drug Acceptors/Donors
bigdf['Norm E-E H-Bonds: DrHA'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm E-E H-Bonds: DrHD'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Drug H-Donors']
bigdf['Norm E-E H-Bonds: TDrAD'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Sum Drug A/D']

#Total Hydrogen Acceptors, Total Hydrogen Donors
bigdf['Norm E-E H-Bonds: THA'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Drug H-Acceptors']
bigdf['Norm E-E H-Bonds: THD'] = bigdf['Total Excipient-Excipient H-Bonds']/ bigdf['Drug H-Donors']

###########################################################################################

#%% Make the percentage N-N Bonds Column

bigdf['Percent N-N Bonds'] =  bigdf['Total N-N H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage N-O Bonds Column

bigdf['Percent N-O Bonds'] =  bigdf['Total N-O H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage O-N Bonds Column

bigdf['Percent O-N Bonds'] =  bigdf['Total O-N H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage O-O Bonds Column

bigdf['Percent O-O Bonds'] =  bigdf['Total O-O H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage F-O Bonds Column

bigdf['Percent F-O Bonds'] =  bigdf['Total F-O H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage F-N Bonds Column

bigdf['Percent F-N Bonds'] =  bigdf['Total F-N H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage O-F Bonds Column

bigdf['Percent O-F Bonds'] =  bigdf['Total O-F H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage N-F Bonds Column

bigdf['Percent N-F Bonds'] =  bigdf['Total N-F H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


###########################################################################################

#%% Make the percentage F-F Bonds Column

bigdf['Percent F-F Bonds'] =  bigdf['Total F-F H-Bonds'] / (bigdf['Total Excipient-Drug H-Bonds'] + bigdf['Total Drug-Excipient H-Bonds'])


############################################################################################

#%% Output DataFrame to .csv
bigdf.to_csv(args.name + '_HBondDF.csv')





