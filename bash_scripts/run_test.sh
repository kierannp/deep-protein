#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH --job-name=test_mpi
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=day-long-std
#SBATCH --exclude=node7,node9,node10,node3
#SBATCH --output=%J.txt

SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep-protein/solvation
cd $SLURM_SUBMIT_DIR

module load gromacs

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
time gmx mdrun -deffnm nvt
time gmx mdrun -deffnm nvt -nb gpu
time gmx mdrun -nb gpu -pme gpu -deffnm nvt 
time srun -n 2 gmx mdrun -deffnm nvt


