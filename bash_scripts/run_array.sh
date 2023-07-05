#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=grid_%J.err
#SBATCH --job-name=des
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=week-long-std
#SBATCH --output=grid_%J.txt
#SBATCH --array=1-4



SLURM_SUBMIT_DIR=/raid6/homes/kierannp/projects/deep-protein
cd $SLURM_SUBMIT_DIR

module load gromacs
module load anaconda/3.9
conda activate deep-protein

# mkdir multi_grid_$WORKER_NUM
# cp bash_scripts/run_grid.sh multi_grid_$WORKER_NUM/grid_$WORKER_NUM.sh
# cd multi_grid_$WORKER_NUM
# # sbatch -J des_$WORKER_NUM grid_$WORKER_NUM.sh --export=WORKER_NUM
# cd ..
python run_grid.py --worker_num $SLURM_ARRAY_TASK_ID

