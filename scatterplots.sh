#!/usr/bin/bash

# remember I took out the dollar sign needed for IPS
# -m be
# -M m.flores@bcbl.eu
# -S /bin/bash

module load python/python3.9

source /bcbl/home/home_g-m/mflores/conda_envs/plotnine_venv/bin/activate

# we are using rys script mprageize_bids.py 
python /bcbl/home/public/MarcoMotion/scripts/Thermal-denoising/scatters_brain.py --source_directory /bcbl/home/public/MarcoMotion/Resting_State/analysis  --subjects sub-001 sub-002 sub-003 sub-004 sub-005 --tasks task-HABLA1200 task-HABLA1700 --methods vanilla nordic hydra tmmpca mppca nordic_magn --overwrite True
