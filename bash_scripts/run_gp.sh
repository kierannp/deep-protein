#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=gp/%J.err
#SBATCH --job-name=des_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=week-long-std
#SBATCH --output=gp/%J.txt



SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep-protein
cd $SLURM_SUBMIT_DIR

module load gromacs
module load anaconda/3.9
conda activate deep-protein

python run_system.py
