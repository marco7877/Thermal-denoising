#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=NORDIC_RS_WOmag
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --mem=80G
#SBATCH --array=1-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt


module load AFNI/24.3.06-foss-2022a
module load MATLAB-IKUR



echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1200" "task-HABLA1700" )
#runs=( "1" "2" )
list_subj=( "001" "002" "003" "004" "005" )

counter=0
for n_sub in "${list_subj[@]}" 
do
	for task in "${tasks[@]}"
	do
#		for run in "${runs[@]}"
#		do
			((counter++))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 2

			fi
		done
	done
#done


method=nordic # nordic, hydra 

subj=sub-${n_sub}
# define folders 
preproc=func_preproc_${method}
repo=/home/mflores/scripts
origin=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}
output=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}_magn

noise=3 # number of noise volumes


save_gfactor=1 # save gfactor map = 1 (saves gfactor and runs nordic); = 2 (saves gfactor and stops)  
# MAIN

# mkdir -p gonna make a chain of dirextories.  
#if [[ ! -e ${output} ]]; then
echo "##################################################"
echo "############# SUBJECT :${subj}      #############"
echo "############# TASK :${task}          #############"
echo "############# RUN :${run}            #############"
echo "##################################################"


#fi


mkdir -p ${output}/gctmp # define temporary folder to save temporary matlab scripts 

# TODO: Change nordic.m to return compressed files them)
list_echoes=$( count -digits 1 1 4 )
echo "******************************************************************************"
for n_echo in ${list_echoes}
do

	temp_nordic=$(mktemp ${output}/gctmp/nordicXXXXXX.m)

	part_mag=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}_dsd

	part_mag_out=${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}N

	part_phase=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-phase_bold_${method}_dsd

	echo "##################################################"
	echo "##################################################"
	echo "##################################################"
	echo "my magnitude volume is: ${part_mag}"

	echo "my output magnitude volume is: ${part_mag_out}"

	echo "my nordic file is in: ${temp_nordic}"

	echo "addpath(genpath('/data/mflores/nordic/'));cd '${repo}';ARG.noise_volume_last=${noise};ARG.temporal_phase=1;ARG.phase_filter_width=10;ARG.save_add_info=1;ARG.magnitude_only=1;ARG.save_gfactor_map=${save_gfactor};ARG.DIROUT='${output}/';NIFTI_NORDIC('${part_mag}.nii.gz','${part_phase}.nii.gz','${part_mag_out}',ARG);" > $temp_nordic

	echo "##################################################"
	echo "##################################################"
	echo "##################################################"
	echo "Thermal denoising with NORDIC:  ID ${subj}"
		matlab -batch "run('$temp_nordic');exit" 

		gzip -f ${output}/${part_mag_out}.nii
done
#rm -rf ${output}/gctmp

echo "******************************************************************************"

volumes=$(3dinfo -nvi ${output}/${part_mag_out}.nii.gz)

echo "3dinfo -nvi ${output}/${part_mag_out}.nii.gz"

vol=$(($volumes-$noise))

echo "Original data ${volumes} volumes, but after removing noise volumes it has ${vol} volumes"

echo "******************************************************************************"


for n_echo in ${list_echoes}
do
	part_mag_thrm=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}N

	part_mag_thrm_trimmed=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}

	if [[ ! -e ${part_mag_thrm_trimmed} ]]; then

		echo "******************************************************************************"
		echo "removing noise volumes to ${part_mag_thrm}.nii"

		3dcalc -a "${part_mag_thrm}.nii.gz[0..${vol}]" -expr 'a' -prefix ${part_mag_thrm_trimmed}.nii.gz -overwrite

		echo "******************************************************************************"
		echo "deleting  temporal ${part_mag_thrm}"
#		rm ${part_mag_thrm}
	fi
done
