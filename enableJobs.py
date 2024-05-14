import os


for x in range(1,131):
    os.system('sbatch AutoJob{}.sh'.format(x))

