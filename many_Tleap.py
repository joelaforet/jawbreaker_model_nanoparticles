import os
import glob


files = glob.glob("ab*")

TLEAP_TEMPLATE = """
source leaprc.gaff
source oldff/leaprc.ff99SB
source leaprc.water.tip3p
loadamberparams Fulvestrant.frcmod
loadamberparams Sorafenib.frcmod
loadamberparams CandesartanCilexetil.frcmod
loadamberparams Taxol.frcmod
loadamberparams CongoRed.frcmod
loadamberparams EvansBlue.frcmod
loadamberparams Lapatinib.frcmod
loadamberparams Glycyrrhizin.frcmod
ZXF = loadmol2 Sorafenib.mol2
ZMH = loadmol2 Fulvestrant.mol2
ZXP = loadmol2 CandesartanCilexetil.mol2
ZPM = loadmol2 Taxol.mol2
ZCR = loadmol2 CongoRed.mol2
ZEB = loadmol2 EvansBlue.mol2
LAP = loadmol2 Lapatinib.mol2
ZGY = loadmol2 Glycyrrhizin.mol2
box = loadPdb %(box)s.pdb
addions box Na+ 0
addions box Cl- 0
saveAmberParm box %(box)s.prmtop %(box)s.inpcrd
quit
"""

for file in files:
    tag = file.split(".")[0]
    tleap_commands = TLEAP_TEMPLATE % dict(box=tag)
    file_handle = open(f'tleap_commands_{tag}.txt', 'w')
    file_handle.writelines(tleap_commands)
    file_handle.close()
    os.system(f"python slurm.py -i 'tleap -f tleap_commands_{tag}.txt' -r -n {tag}") 

print("Job's Done!")

