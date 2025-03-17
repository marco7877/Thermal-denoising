#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=align_mppca_RS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:40:00
#SBATCH --array=1-10
#SBATCH --mem=16G
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI 
module load FSL

tasks=( "task-HABLA1700" "task-HABLA1200")
list_subj=( "sub-001" "sub-002" "sub-003" "sub-004" "sub-005" )
#tasks=( "task-HABLA1700" )
#list_subj=( "sub-003" )

counter=0
for subj in "${list_subj[@]}" 
do
	for tsk in "${tasks[@]}"
	do
			((counter++))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 2

			fi
	done
done

echo "##################################################"
echo "############# SUBJECT :${subj}      #############"
echo "############# TASK :${tsk}          #############"
echo "##################################################"


method=mppca # vanilla, nordic, hydra 
#n_subj=${n_sub}
task=${tsk}
preproc=func_preproc_${method}

# define folders
#subj=sub-00${n_subj}
origin=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}
output=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}
input=/scratch/mflores/Resting_State/${subj}/ses-1/func

nechoes=4


# DSD extention for the volumes without thermal noise (nordic, hydra) or vainilla (ME)

ext=.nii.gz
# MAIN

#Check if output durectory exists

if [[ ! -e ${output} ]]; then

	mkdir ${output}

fi

# Reference original volumnes
part_mag=${origin}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}${ext}
echo "My reference part_mag is: ${part_mag}"
sbref=${input}/${subj}_ses-1_${task}_echo-1_part-mag_sbref${ext}
# define name for symbolic link 
sbref_ln=${output}/${subj}_ses-1_${task}_echo-1_part-mag_sbref${ext}

if [[ ! -e ${sbref_ln} ]]; then
	ln -s ${sbref} ${sbref_ln}
fi

#TODO. Add indentation with | to split loong commands
#  Make variables more clean. let only the extension after variable, i.e., ${output_mcf}_al-nii.gz
3dvolreg -overwrite -Fourier -base ${sbref_ln} \
	-1Dfile  ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
	-1Dmatrix_save ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D \
	-prefix ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf${ext} \
	${part_mag}

#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
#		-demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D \
#		-derivative -demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_deriv1.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
# 		-derivative -collapse_cols euclidean_norm -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_enorm.1D -overwrite

list_echoes=( "1" "2" "3" "4" )

for n_echo in "${list_echoes[@]}"
do

	part_mag=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${ext}

	echo "3dAllineate -overwrite -base ${sbref_ln} -final cubic -1Dmatrix_apply ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}_mcf_al${ext} ${part_mag}"

	3dAllineate -overwrite -base ${sbref_ln} -final cubic \
		-1Dmatrix_apply ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D \
		-prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}_mcf_al${ext} \
		${part_mag} 

	
done
