#!/bin/bash
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jrl78@duke.edu
#SBATCH -e slurm.err
#SBATCH --job-name=10ns_DNA.job
#SBATCH --mem=20G
#SBATCH -p scavenger-gpu --gres=gpu:1
#SBATCH --exclusive

/bin/csh SIMULATE
