import os
import os.path
import mdtraj as md
import pandas as pd
import openpyxl as pxl

n_frames = 200
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        Hdonors = {}
        Hacceptors = {}
        
        for i in range(n_frames):
            traj = md.load_frame(filename, i)
            hbonds = md.baker_hubbard(traj)
            for hbond in hbonds:
                r1 = str(traj.topology.atom(hbond[0]))
                r1 = r1[0:3] + '-' + (r1[-1] if r1[-1].isnumeric() is False else r1[-2:]) 
                r2 = str(traj.topology.atom(hbond[2]))
                r2 = r2[0:3] + '-' + (r2[-1] if r2[-1].isnumeric() is False else r2[-2:])
                if ('ZPE' in r1) is not ('ZPE' in r2):
                    Hdonors[r1] = 1 if r1 not in Hdonors else Hdonors[r1] + 1
                    Hacceptors[r2] = 1 if r2 not in Hacceptors else Hacceptors[r2] + 1
                    
        df = pd.DataFrame({'Hdonors': Hdonors, 'Hacceptors': Hacceptors})
        
        if not os.path.isfile('Hbonds_count.xlsx'):
            df.to_excel('Hbonds_count.xlsx', '%s' % excipient_name, index = True)

        else:
            excel_book = pxl.load_workbook('Hbonds_count.xlsx')
            with pd.ExcelWriter('Hbonds_count.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
                
            

        
