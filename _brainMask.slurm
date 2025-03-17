#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=BrainMask_RS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=30G
#SBATCH --array=1-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt


module load Anaconda3

module load AFNI
sub=${SLURM_ARRAY_TASK_ID}
python ~/scripts/brainMask.py --bids_dir /scratch/mflores/Resting_State --filt_pattern sub-00${sub} 
