#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=grid_%J.err
#SBATCH --job-name=des_$1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=week-long-std
#SBATCH --output=grid_%J.txt



SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep-protein
cd $SLURM_SUBMIT_DIR

module load gromacs
module load anaconda/3.9
conda activate deep-protein

echo "The number is $WORKER_NUM"

python run_grid.py --worker_num $WORKER_NUM
