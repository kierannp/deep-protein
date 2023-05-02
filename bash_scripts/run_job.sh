#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH --job-name=des_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=day-long-std
#SBATCH --exclude=node7,node9,node10,node3
#SBATCH --output=%J.txt

SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep-protein/solvation
cd $SLURM_SUBMIT_DIR

module load gromacs

gmx grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -v -deffnm md_0_1 -nt 4
