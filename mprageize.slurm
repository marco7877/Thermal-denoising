#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=MPRAGE_RS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI
module load Anaconda3


# we are using rys script mprageize_bids.py 

python ~/scripts/mprageize_bids.py --bids_dir /scratch/mflores/Resting_State  --mprageize_dir /scratch/mflores/myToolboxes/3dMPRAGEise --overwrite True
