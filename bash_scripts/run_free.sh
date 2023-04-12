#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=free_%J.err
#SBATCH --output=free_%J.out
#SBATCH --job-name=free_energy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=day-long-std
#SBATCH --exclude=node7,node9,node10,node3
#SBATCH --output=free_%J.txt

SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep_protein/free-energy
cd $SLURM_SUBMIT_DIR

module load gromacs

bash job.sh
