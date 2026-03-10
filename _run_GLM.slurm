#!/bin/bash 
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=GLMtask_predict_
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --begin=now
#SBATCH --mem=800G
#SBATCH --array=1-5
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/MEVOTLOC/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/MEVOTLOC/logs/%x_%a_%j.txt

module load Miniforge3
source activate /scratch/mflores/nilearn_emerald

methods=("vanilla" "mppca" "nordic" "tmppca" "hydra")

#methods=("vanilla")
counter=0
for method in "${methods[@]}"
do
	((counter++))
	echo "${counter}"
	echo "${SLURM_ARRAY_TASK_ID}"
	if [[ "${counter}" -eq "${SLURM_ARRAY_TASK_ID}" ]]; then
		echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
		break 1
	fi
done
echo "##################################################"
echo "############# METHOD :${method}            #############"
echo "##################################################"

cd /home/mflores/Thermal-denoising
func=func_preproc_${method}

dir1=/data/mflores/MEVOTLOC/sub-001/ses-01/${func}/tedana_task-fLoc_run
dir2=/data/mflores/MEVOTLOC/sub-001/ses-02/${func}/tedana_task-fLoc_run
dir3=/data/mflores/MEVOTLOC/sub-001/ses-03/${func}/tedana_task-fLoc_run
dir4=/data/mflores/MEVOTLOC/sub-001/ses-04/${func}/tedana_task-fLoc_run

mkdir /scratch/mflores/MEVOTLOC/${method}_task


python _GLM.py --nii_files ${dir1}-1_meicabold/desc-optcom_bold.nii.gz ${dir1}-2_meicabold/desc-optcom_bold.nii.gz ${dir1}-4_meicabold/desc-optcom_bold.nii.gz ${dir1}-6_meicabold/desc-optcom_bold.nii.gz ${dir1}-7_meicabold/desc-optcom_bold.nii.gz ${dir2}-6_meicabold/desc-optcom_bold.nii.gz ${dir2}-7_meicabold/desc-optcom_bold.nii.gz ${dir2}-8_meicabold/desc-optcom_bold.nii.gz ${dir2}-9_meicabold/desc-optcom_bold.nii.gz ${dir3}-1_meicabold/desc-optcom_bold.nii.gz ${dir3}-3_meicabold/desc-optcom_bold.nii.gz ${dir3}-4_meicabold/desc-optcom_bold.nii.gz ${dir3}-5_meicabold/desc-optcom_bold.nii.gz ${dir3}-6_meicabold/desc-optcom_bold.nii.gz ${dir3}-7_meicabold/desc-optcom_bold.nii.gz ${dir3}-10_meicabold/desc-optcom_bold.nii.gz ${dir4}-1_meicabold/desc-optcom_bold.nii.gz ${dir4}-4_meicabold/desc-optcom_bold.nii.gz ${dir4}-6_meicabold/desc-optcom_bold.nii.gz ${dir4}-7_meicabold/desc-optcom_bold.nii.gz ${dir4}-8_meicabold/desc-optcom_bold.nii.gz --events_files /data/mflores/MEVOTLOC/designMatrix_01.csv /data/mflores/MEVOTLOC/designMatrix_02.csv /data/mflores/MEVOTLOC/designMatrix_04.csv /data/mflores/MEVOTLOC/designMatrix_06.csv /data/mflores/MEVOTLOC/designMatrix_07.csv /data/mflores/MEVOTLOC/designMatrix_06.csv /data/mflores/MEVOTLOC/designMatrix_07.csv /data/mflores/MEVOTLOC/designMatrix_08.csv /data/mflores/MEVOTLOC/designMatrix_09.csv /data/mflores/MEVOTLOC/designMatrix_01.csv /data/mflores/MEVOTLOC/designMatrix_03.csv /data/mflores/MEVOTLOC/designMatrix_04.csv /data/mflores/MEVOTLOC/designMatrix_05.csv /data/mflores/MEVOTLOC/designMatrix_06.csv /data/mflores/MEVOTLOC/designMatrix_07.csv /data/mflores/MEVOTLOC/designMatrix_10.csv /data/mflores/MEVOTLOC/designMatrix_01.csv /data/mflores/MEVOTLOC/designMatrix_04.csv /data/mflores/MEVOTLOC/designMatrix_06.csv /data/mflores/MEVOTLOC/designMatrix_07.csv /data/mflores/MEVOTLOC/designMatrix_08.csv --mask /scratch/mflores/MEVOTLOC/masks/sub-001_ses-01_task-fLoc_run-1_echo-1_part-mag_gm-mask.nii.gz --n_task_regressors 6 --n_splits 100 --test_size 0.3 --n_permutations 0 --random_state 42 --debug --output_prefix /scratch/mflores/MEVOTLOC/${method}_task
