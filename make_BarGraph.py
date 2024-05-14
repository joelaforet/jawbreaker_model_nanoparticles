"""
@author: Joe Laforet
Created on Thu May 27 16:53:12 2021


I understand and have adhered to all the tenets of the Duke Community Standard
in creating this code.
Signed: [jrl78]

"""
#%% Import modules
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

#%% Initialize Argparse

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument('-i', '--infile', help = 'Input .csv of HBond data for a simulation set.', required = True)
parser.add_argument('-n', '--name', help = "Input name for output files.", required = True)
parser.add_argument('-d', '--directory', help = "Input directory where data is stored.", required = True)
parser.add_argument('-f', '--foldername', help = 'Input name for folder where bar graphs will be generated in.',  default = 'BarGraphs')

args = parser.parse_args()

#%%

os.chdir(args.directory)
curDir = os.getcwd()

path = args.foldername

if not os.path.exists(path):
    os.makedirs(path)


f = args.infile
indf = pd.read_csv(args.infile)
indf.drop(columns = ['Unnamed: 0'], inplace = True)

os.chdir(path)

#%% Make the bar graphs

#######################################################################################

def getName(row):
    return row.split('_')[1]

#%% Make bar graph of total H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total HBonds formed')
#indf['Temporary'] = indf['Tag'].apply(getName)
#print(indf['Temporary'])
indf.sort_values(by = 'Total H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total H-Bonds'])

fig.tight_layout()
fig.set_size_inches(8,6, forward = True)
fig.savefig(('{}_TotalHBonds.pdf').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Excipient Acceptors/Donors)
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: TExAD')
indf.sort_values(by = 'Norm Tot H-Bonds: TExAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: TExAD'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_TotalHBonds_TExAD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Drug Acceptors/Donors)
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig1.add_subplot(1,1,1)
ax2.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: TDrAD')
indf.sort_values(by = 'Norm Tot H-Bonds: TDrAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: TDrAD'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_TotalHBonds_TDrAD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Hydrogen Acceptors in both species)
fig3 = plt.figure(num = 1, clear = True)
ax3 = fig1.add_subplot(1,1,1)
ax3.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: THA')
indf.sort_values(by = 'Norm Tot H-Bonds: THA', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: THA'])

fig3.tight_layout()
fig3.set_size_inches(6,4, forward = True)
fig3.savefig(('{}_TotalHBonds_THA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Hydrogen Donors in both species)
fig4 = plt.figure(num = 1, clear = True)
ax4 = fig1.add_subplot(1,1,1)
ax4.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: THD')
indf.sort_values(by = 'Norm Tot H-Bonds: THD', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: THD'])

fig4.tight_layout()
fig4.set_size_inches(6,4, forward = True)
fig4.savefig(('{}_TotalHBonds_THD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Excipient Hydrogen Acceptors)
fig5 = plt.figure(num = 1, clear = True)
ax5 = fig1.add_subplot(1,1,1)
ax5.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: ExHA')
indf.sort_values(by = 'Norm Tot H-Bonds: ExHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: ExHA'])

fig5.tight_layout()
fig5.set_size_inches(6,4, forward = True)
fig5.savefig(('{}_TotalHBonds_ExHA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Excipient Hydrogen Donors)
fig6 = plt.figure(num = 1, clear = True)
ax6 = fig1.add_subplot(1,1,1)
ax6.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: ExHD')
indf.sort_values(by = 'Norm Tot H-Bonds: ExHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: ExHD'])

fig6.tight_layout()
fig6.set_size_inches(6,4, forward = True)
fig6.savefig(('{}_TotalHBonds_ExHD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Drug Hydrogen Acceptors)
fig7 = plt.figure(num = 1, clear = True)
ax7 = fig1.add_subplot(1,1,1)
ax7.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: DrHA')
indf.sort_values(by = 'Norm Tot H-Bonds: DrHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: DrHA'])

fig7.tight_layout()
fig7.set_size_inches(6,4, forward = True)
fig7.savefig(('{}_TotalHBonds_DrHA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Drug Hydrogen Donors)
fig8 = plt.figure(num = 1, clear = True)
ax8 = fig1.add_subplot(1,1,1)
ax8.set(ylabel = 'Drug/Excipient Pair', xlabel = 'HBonds formed', title = 'Total H-Bonds Normalized: DrHD')
indf.sort_values(by = 'Norm Tot H-Bonds: DrHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm Tot H-Bonds: DrHD'])

fig8.tight_layout()
fig8.set_size_inches(6,4, forward = True)
fig8.savefig(('{}_TotalHBonds_DrHD.png').format(args.name))

#######################################################################################

#%% Make bar graph of total D-D H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D HBonds formed')
indf.sort_values(by = 'Total Drug-Drug H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total Drug-Drug H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_ddHBonds.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Excipient Acceptors/Donors)
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: TExAD')
indf.sort_values(by = 'Norm D-D H-Bonds: TExAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: TExAD'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_ddHBonds_TExAD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Drug Acceptors/Donors)
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig1.add_subplot(1,1,1)
ax2.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: TDrAD')
indf.sort_values(by = 'Norm D-D H-Bonds: TDrAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: TDrAD'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_ddHBonds_TDrAD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Hydrogen Acceptors in both species)
fig3 = plt.figure(num = 1, clear = True)
ax3 = fig1.add_subplot(1,1,1)
ax3.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: THA')
indf.sort_values(by = 'Norm D-D H-Bonds: THA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: THA'])

fig3.tight_layout()
fig3.set_size_inches(6,4, forward = True)
fig3.savefig(('{}_ddHBonds_THA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Sum of Hydrogen Donors in both species)
fig4 = plt.figure(num = 1, clear = True)
ax4 = fig1.add_subplot(1,1,1)
ax4.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: THD')
indf.sort_values(by = 'Norm D-D H-Bonds: THD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: THD'])

fig4.tight_layout()
fig4.set_size_inches(6,4, forward = True)
fig4.savefig(('{}_ddHBonds_THD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Excipient Hydrogen Acceptors)
fig5 = plt.figure(num = 1, clear = True)
ax5 = fig1.add_subplot(1,1,1)
ax5.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: ExHA')
indf.sort_values(by = 'Norm D-D H-Bonds: ExHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: ExHA'])

fig5.tight_layout()
fig5.set_size_inches(6,4, forward = True)
fig5.savefig(('{}_ddHBonds_ExHA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Excipient Hydrogen Donors)
fig6 = plt.figure(num = 1, clear = True)
ax6 = fig1.add_subplot(1,1,1)
ax6.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: ExHD')
indf.sort_values(by = 'Norm D-D H-Bonds: ExHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: ExHD'])

fig6.tight_layout()
fig6.set_size_inches(6,4, forward = True)
fig6.savefig(('{}_ddHBonds_ExHD.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Drug Hydrogen Acceptors)
fig7 = plt.figure(num = 1, clear = True)
ax7 = fig1.add_subplot(1,1,1)
ax7.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: DrHA')
indf.sort_values(by = 'Norm D-D H-Bonds: DrHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: DrHA'])

fig7.tight_layout()
fig7.set_size_inches(6,4, forward = True)
fig7.savefig(('{}_ddHBonds_DrHA.png').format(args.name))

#%% Make bar graph of total H-Bonds (Normalized for Drug Hydrogen Donors)
fig8 = plt.figure(num = 1, clear = True)
ax8 = fig1.add_subplot(1,1,1)
ax8.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-D HBonds formed', title = 'Total D-D H-Bonds Normalized: DrHD')
indf.sort_values(by = 'Norm D-D H-Bonds: DrHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-D H-Bonds: DrHD'])

fig8.tight_layout()
fig8.set_size_inches(6,4, forward = True)
fig8.savefig(('{}_ddHBonds_DrHD.png').format(args.name))

#######################################################################################

#%% Make bar graph of total D-E H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E HBonds formed')
indf.sort_values(by = 'Total Drug-Excipient H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total Drug-Excipient H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_deHBonds.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Sum of Excipient Acceptors/Donors)
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: TExAD')
indf.sort_values(by = 'Norm D-E H-Bonds: TExAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: TExAD'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_deHBonds_TExAD.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Sum of Drug Acceptors/Donors)
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig1.add_subplot(1,1,1)
ax2.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: TDrAD')
indf.sort_values(by = 'Norm D-E H-Bonds: TDrAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: TDrAD'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_deHBonds_TDrAD.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Sum of Hydrogen Acceptors in both species)
fig3 = plt.figure(num = 1, clear = True)
ax3 = fig1.add_subplot(1,1,1)
ax3.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: THA')
indf.sort_values(by = 'Norm D-E H-Bonds: THA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: THA'])

fig3.tight_layout()
fig3.set_size_inches(6,4, forward = True)
fig3.savefig(('{}_deHBonds_THA.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Sum of Hydrogen Donors in both species)
fig4 = plt.figure(num = 1, clear = True)
ax4 = fig1.add_subplot(1,1,1)
ax4.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: THD')
indf.sort_values(by = 'Norm D-E H-Bonds: THD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: THD'])

fig4.tight_layout()
fig4.set_size_inches(6,4, forward = True)
fig4.savefig(('{}_deHBonds_THD.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Excipient Hydrogen Acceptors)
fig5 = plt.figure(num = 1, clear = True)
ax5 = fig1.add_subplot(1,1,1)
ax5.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: ExHA')
indf.sort_values(by = 'Norm D-E H-Bonds: ExHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: ExHA'])

fig5.tight_layout()
fig5.set_size_inches(6,4, forward = True)
fig5.savefig(('{}_deHBonds_ExHA.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Excipient Hydrogen Donors)
fig6 = plt.figure(num = 1, clear = True)
ax6 = fig1.add_subplot(1,1,1)
ax6.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: ExHD')
indf.sort_values(by = 'Norm D-E H-Bonds: ExHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: ExHD'])

fig6.tight_layout()
fig6.set_size_inches(6,4, forward = True)
fig6.savefig(('{}_deHBonds_ExHD.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Drug Hydrogen Acceptors)
fig7 = plt.figure(num = 1, clear = True)
ax7 = fig1.add_subplot(1,1,1)
ax7.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: DrHA')
indf.sort_values(by = 'Norm D-E H-Bonds: DrHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: DrHA'])

fig7.tight_layout()
fig7.set_size_inches(6,4, forward = True)
fig7.savefig(('{}_deHBonds_DrHA.png').format(args.name))

#%% Make bar graph of total D-E H-Bonds (Normalized for Drug Hydrogen Donors)
fig8 = plt.figure(num = 1, clear = True)
ax8 = fig1.add_subplot(1,1,1)
ax8.set(ylabel = 'Drug/Excipient Pair', xlabel = 'D-E HBonds formed', title = 'Total D-E H-Bonds Normalized: DrHD')
indf.sort_values(by = 'Norm D-E H-Bonds: DrHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm D-E H-Bonds: DrHD'])

fig8.tight_layout()
fig8.set_size_inches(6,4, forward = True)
fig8.savefig(('{}_deHBonds_DrHD.png').format(args.name))

#######################################################################################

#%% Make bar graph of total E-D H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D HBonds formed')
indf.sort_values(by = 'Total Excipient-Drug H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total Excipient-Drug H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_edHBonds.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Sum of Excipient Acceptors/Donors)
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: TExAD')
indf.sort_values(by = 'Norm E-D H-Bonds: TExAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: TExAD'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_edHBonds_TExAD.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Sum of Drug Acceptors/Donors)
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig1.add_subplot(1,1,1)
ax2.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: TDrAD')
indf.sort_values(by = 'Norm E-D H-Bonds: TDrAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: TDrAD'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_edHBonds_TDrAD.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Sum of Hydrogen Acceptors in both species)
fig3 = plt.figure(num = 1, clear = True)
ax3 = fig1.add_subplot(1,1,1)
ax3.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: THA')
indf.sort_values(by = 'Norm E-D H-Bonds: THA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: THA'])

fig3.tight_layout()
fig3.set_size_inches(6,4, forward = True)
fig3.savefig(('{}_edHBonds_THA.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Sum of Hydrogen Donors in both species)
fig4 = plt.figure(num = 1, clear = True)
ax4 = fig1.add_subplot(1,1,1)
ax4.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: THD')
indf.sort_values(by = 'Norm E-D H-Bonds: THD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: THD'])

fig4.tight_layout()
fig4.set_size_inches(6,4, forward = True)
fig4.savefig(('{}_edHBonds_THD.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Excipient Hydrogen Acceptors)
fig5 = plt.figure(num = 1, clear = True)
ax5 = fig1.add_subplot(1,1,1)
ax5.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: ExHA')
indf.sort_values(by = 'Norm E-D H-Bonds: ExHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: ExHA'])

fig5.tight_layout()
fig5.set_size_inches(6,4, forward = True)
fig5.savefig(('{}_edHBonds_ExHA.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Excipient Hydrogen Donors)
fig6 = plt.figure(num = 1, clear = True)
ax6 = fig1.add_subplot(1,1,1)
ax6.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: ExHD')
indf.sort_values(by = 'Norm E-D H-Bonds: ExHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: ExHD'])

fig6.tight_layout()
fig6.set_size_inches(6,4, forward = True)
fig6.savefig(('{}_edHBonds_ExHD.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Drug Hydrogen Acceptors)
fig7 = plt.figure(num = 1, clear = True)
ax7 = fig1.add_subplot(1,1,1)
ax7.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: DrHA')
indf.sort_values(by = 'Norm E-D H-Bonds: DrHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: DrHA'])

fig7.tight_layout()
fig7.set_size_inches(6,4, forward = True)
fig7.savefig(('{}_edHBonds_DrHA.png').format(args.name))

#%% Make bar graph of total E-D H-Bonds (Normalized for Drug Hydrogen Donors)
fig8 = plt.figure(num = 1, clear = True)
ax8 = fig1.add_subplot(1,1,1)
ax8.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-D HBonds formed', title = 'Total E-D H-Bonds Normalized: DrHD')
indf.sort_values(by = 'Norm E-D H-Bonds: DrHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-D H-Bonds: DrHD'])

fig8.tight_layout()
fig8.set_size_inches(6,4, forward = True)
fig8.savefig(('{}_edHBonds_DrHD.png').format(args.name))

#######################################################################################

#%% Make bar graph of total E-E H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E HBonds formed')
indf.sort_values(by = 'Total Excipient-Excipient H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total Excipient-Excipient H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_eeHBonds.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Sum of Excipient Acceptors/Donors)
fig1 = plt.figure(num = 1, clear = True)
ax1 = fig1.add_subplot(1,1,1)
ax1.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: TExAD')
indf.sort_values(by = 'Norm E-E H-Bonds: TExAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: TExAD'])

fig1.tight_layout()
fig1.set_size_inches(6,4, forward = True)
fig1.savefig(('{}_eeHBonds_TExAD.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Sum of Drug Acceptors/Donors)
fig2 = plt.figure(num = 1, clear = True)
ax2 = fig1.add_subplot(1,1,1)
ax2.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: TDrAD')
indf.sort_values(by = 'Norm E-E H-Bonds: TDrAD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: TDrAD'])

fig2.tight_layout()
fig2.set_size_inches(6,4, forward = True)
fig2.savefig(('{}_eeHBonds_TDrAD.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Sum of Hydrogen Acceptors in both species)
fig3 = plt.figure(num = 1, clear = True)
ax3 = fig1.add_subplot(1,1,1)
ax3.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: THA')
indf.sort_values(by = 'Norm E-E H-Bonds: THA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: THA'])

fig3.tight_layout()
fig3.set_size_inches(6,4, forward = True)
fig3.savefig(('{}_eeHBonds_THA.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Sum of Hydrogen Donors in both species)
fig4 = plt.figure(num = 1, clear = True)
ax4 = fig1.add_subplot(1,1,1)
ax4.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: THD')
indf.sort_values(by = 'Norm E-E H-Bonds: THD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: THD'])

fig4.tight_layout()
fig4.set_size_inches(6,4, forward = True)
fig4.savefig(('{}_eeHBonds_THD.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Excipient Hydrogen Acceptors)
fig5 = plt.figure(num = 1, clear = True)
ax5 = fig1.add_subplot(1,1,1)
ax5.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: ExHA')
indf.sort_values(by = 'Norm E-E H-Bonds: ExHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: ExHA'])

fig5.tight_layout()
fig5.set_size_inches(6,4, forward = True)
fig5.savefig(('{}_eeHBonds_ExHA.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Excipient Hydrogen Donors)
fig6 = plt.figure(num = 1, clear = True)
ax6 = fig1.add_subplot(1,1,1)
ax6.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: ExHD')
indf.sort_values(by = 'Norm E-E H-Bonds: ExHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: ExHD'])

fig6.tight_layout()
fig6.set_size_inches(6,4, forward = True)
fig6.savefig(('{}_eeHBonds_ExHD.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Drug Hydrogen Acceptors)
fig7 = plt.figure(num = 1, clear = True)
ax7 = fig1.add_subplot(1,1,1)
ax7.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: DrHA')
indf.sort_values(by = 'Norm E-E H-Bonds: DrHA', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: DrHA'])

fig7.tight_layout()
fig7.set_size_inches(6,4, forward = True)
fig7.savefig(('{}_eeHBonds_DrHA.png').format(args.name))

#%% Make bar graph of total E-E H-Bonds (Normalized for Drug Hydrogen Donors)
fig8 = plt.figure(num = 1, clear = True)
ax8 = fig1.add_subplot(1,1,1)
ax8.set(ylabel = 'Drug/Excipient Pair', xlabel = 'E-E HBonds formed', title = 'Total E-E H-Bonds Normalized: DrHD')
indf.sort_values(by = 'Norm E-E H-Bonds: DrHD', inplace = True)
plt.barh(indf['Tag'], indf['Norm E-E H-Bonds: DrHD'])

fig8.tight_layout()
fig8.set_size_inches(6,4, forward = True)
fig8.savefig(('{}_eeHBonds_DrHD.png').format(args.name))

#######################################################################################

#%% Make bar graph of total N-N H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'N-N HBonds formed', title = 'Total N-N HBonds formed')
indf.sort_values(by = 'Total N-N H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total N-N H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_NNHBonds.png').format(args.name))

#%% Make bar graph of total N-O H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'N-O HBonds formed', title = 'Total N-O HBonds formed')
indf.sort_values(by = 'Total N-O H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total N-O H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_NOHBonds.png').format(args.name))

#%% Make bar graph of total O-N H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'O-N HBonds formed', title = 'Total O-N HBonds formed')
indf.sort_values(by = 'Total O-N H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total O-N H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_ONHBonds.png').format(args.name))

#%% Make bar graph of total O-O H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'O-O HBonds formed', title = 'Total O-O HBonds formed')
indf.sort_values(by = 'Total O-O H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total O-O H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_OOHBonds.png').format(args.name))

#%% Make bar graph of total F-O H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'F-O HBonds formed', title = 'Total F-O HBonds formed')
indf.sort_values(by = 'Total F-O H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total F-O H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_FOHBonds.png').format(args.name))

#%% Make bar graph of total F-N H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'F-N HBonds formed', title = 'Total F-N HBonds formed')
indf.sort_values(by = 'Total F-N H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total F-N H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_FNHBonds.png').format(args.name))

#%% Make bar graph of total O-F H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'O-F HBonds formed', title = 'Total O-F HBonds formed')
indf.sort_values(by = 'Total O-F H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total O-F H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_OFHBonds.png').format(args.name))

#%% Make bar graph of total N-F H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'N-F HBonds formed', title = 'Total N-F HBonds formed')
indf.sort_values(by = 'Total N-F H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total N-F H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_NFHBonds.png').format(args.name))

#%% Make bar graph of total F-F H-Bonds (Unnormalized)
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1,1,1)
ax.set(ylabel = 'Drug/Excipient Pair', xlabel = 'F-F HBonds formed', title = 'Total F-F HBonds formed')
indf.sort_values(by = 'Total F-F H-Bonds', inplace = True)
plt.barh(indf['Tag'], indf['Total F-F H-Bonds'])

fig.tight_layout()
fig.set_size_inches(6,4, forward = True)
fig.savefig(('{}_FFHBonds.png').format(args.name))


#######################################################################################
