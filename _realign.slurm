#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=align_nordic_residuals
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


tasks=("task-HABLA1700" "task-HABLA1200")
list_subj=("sub-001" "sub-002" "sub-003" "sub-004" "sub-005")
#methods=( "tmmpca" "hydra" "mppca" "vanilla" )
methods=( "nordic" )
# In case of halves, rememberthat extension == ?.nii.gz
ext=.nii.gz
echo "Your files will have the extention: ${ext}"
counter=0
nechoes=4
residuals="True"
trans_matrix="False"

echo My task ID is"${SLURM_ARRAY_TASK_ID}"

for subj in "${list_subj[@]}"; do
	for tsk in "${tasks[@]}"; do
		for method in "${methods[@]}"; do
			((counter++))
			echo "${counter}"
			echo "${SLURM_ARRAY_TASK_ID}"
			if [[ "${counter}" -eq "${SLURM_ARRAY_TASK_ID}" ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 3

			fi
		done
	done
done

echo "##################################################"
echo "############# SUBJECT :${subj}      #############"
echo "############# TASK :${tsk}          #############"
echo "############# METHOD :${method}          #############"
echo "############# RESIDUALS == ${residuals}          #############"
echo "############# calculating matrix == :${trans_matrix}          #############"
echo "##################################################"


task=${tsk}
preproc=func_preproc_${method}

# define folders
origin=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}
output=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}
matrix_dir=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_vanilla
input=/data/HABLA/Resting_State/${subj}/ses-1/func



if [[ ! -e ${output} ]]; then

	mkdir ${output}

fi

# Reference original volumnes
part_mag=${origin}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}${ext}
echo "My reference part_mag is: ${part_mag}"
sbref=${input}/${subj}_ses-1_${task}_echo-1_part-mag_sbref.nii.gz
# define name for symbolic link 
sbref_ln=${output}/${subj}_ses-1_${task}_echo-1_part-mag_sbref.nii.gz

if [[ ! -e ${sbref_ln} ]]; then
	cp -s ${sbref} ${sbref_ln}
fi

matrix=${matrix_dir}/${subj}_ses-1_${task}_echo-1_part-mag_bold_vanilla_mcf.aff12.1D
if [[ ${trans_matrix} == "True" ]];then
	3dvolreg -overwrite -Fourier -base ${sbref_ln} \
		-1Dfile  ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
		-1Dmatrix_save ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D \
		-prefix ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.nii.gz \
		${part_mag}
	matrix=${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D
fi
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
#		-demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D \
#		-derivative -demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_deriv1.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
# 		-derivative -collapse_cols euclidean_norm -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_enorm.1D -overwrite


list_echoes=( "1" "2" "3" "4" )

for n_echo in "${list_echoes[@]}"
do
	if [[ ${residuals} == "True" ]];then
		original=${matrix_dir}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_vanilla${ext}
		processed=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${ext}
		part_mag=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_residuals_${method}${ext}
		aligned=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_residuals_${method}_mcf_al${ext}
		3dcalc -a ${original} -b ${processed} -expr "(a-b)" -prefix ${part_mag} -overwrite 
	else
		part_mag=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${ext}
		aligned=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}_mcf_al${ext}
	fi
	echo "3dAllineate -overwrite -base ${sbref_ln} -final cubic -1Dmatrix_apply ${matrix} -prefix ${aligned} ${part_mag}"

	3dAllineate -overwrite -base ${sbref_ln} -final cubic \
		-1Dmatrix_apply ${matrix} \
		-prefix ${aligned} \
		${part_mag} 

done
