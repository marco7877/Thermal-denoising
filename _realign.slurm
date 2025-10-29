#!/bin/bash -xv
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=align_Rest_HighRes_complete_run-2_vanilla
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --array=2
#SBATCH --mem=32G
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt

PS4='${LINENO}: '

module load AFNI 


tasks=("task-REST" )
list_subj=( "sub-001")
#methods=( "mppca" "nordic" "hydra" "tmppca" )
methods=( "vanilla" )
#methods=("tmppca")
runs=("1" "2")
#runs=( "1")
multiple_runs="True"
# In case of halves, rememberthat extension == ?.nii.gz
counter=0
nechoes=4
residuals="True"
trans_matrix="True"
halves="False"

echo My task ID is"${SLURM_ARRAY_TASK_ID}"

for subj in "${list_subj[@]}"; do
	for tsk in "${tasks[@]}"; do
		for method in "${methods[@]}"; do
#			for run in "${runs[@]}"; do
				((counter++))
				echo "${counter}"
				echo "${SLURM_ARRAY_TASK_ID}"
				if [[ "${counter}" -eq "${SLURM_ARRAY_TASK_ID}" ]]; then
					echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
					break 3

				fi
#			done
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

if [[ ${multiple_runs} == "True" ]];then
task=${tsk}_run-${run}
sbref=${input}/${subj}_ses-1_${tsk}_run-1_echo-1_part-mag_sbref.nii.gz
else
task=${tsk}
sbref=${input}/${subj}_ses-1_${tsk}_echo-1_part-mag_sbref.nii.gz
fi
preproc=func_preproc_${method}

# define folders
origin=/data/mflores/Rest_HighRes/${subj}/ses-1/${preproc}
output=/scratch/mflores/Rest_HighRes/${subj}/ses-1/${preproc}
matrix_dir=/scratch/mflores/Rest_HighRes/${subj}/ses-1/func_preproc_vanilla
input=/data/mflores/Rest_HighRes/${subj}/ses-1/func


mkdir -p ${output}

# Reference original volumnes
echo "My reference part_mag is: ${part_mag}"
# define name for symbolic link 
sbref_ln=${output}/${subj}_ses-1_${task}_echo-1_part-mag_sbref.nii.gz

if [[ ! -e ${sbref_ln} ]]; then
	cp -s ${sbref} ${sbref_ln}
fi

if [[ $halves == "True" ]]; then
	extentions=( "" "1" "2" )
else
	extentions=""
fi
if [[ ${trans_matrix} == "True" ]];then
	for ext in "${extentions[@]}"; do
		part_mag=${origin}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}${ext}.nii.gz
		3dvolreg -overwrite -Fourier -base ${sbref} \
			-1Dfile  ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf${ext}.1D \
			-1Dmatrix_save ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf${ext}.aff12.1D \
			-prefix ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf${ext}.nii.gz \
			${part_mag}
		#	matrix=${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.aff12.1D
	done
fi
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
	#		-demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_demean.1D \
	#		-derivative -demean -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_deriv1.1D -overwrite
#	1d_tool.py -infile ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf.1D \
	# 		-derivative -collapse_cols euclidean_norm -write ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_mcf_enorm.1D -overwrite

list_echoes=($(seq 1 ${nechoes}))

for ext in "${extentions[@]}"; do

	matrix=${matrix_dir}/${subj}_ses-1_${task}_echo-1_part-mag_bold_vanilla_mcf${ext}.aff12.1D
	echo "Your files will have the extention: ${ext}"
	for n_echo in "${list_echoes[@]}"
	do
		if [[ ${residuals} == "True" ]];then
			original=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_vanilla${ext}.nii.gz
			processed=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${ext}.nii.gz
			part_mag=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_residuals_${method}${ext}.nii.gz
			aligned=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_residuals_${method}_mcf_al${ext}.nii.gz
			3dcalc -a ${original} -b ${processed} -expr "(a-b)" -prefix ${part_mag} -overwrite 
			echo "3dAllineate -overwrite -base ${sbref} -final cubic -1Dmatrix_apply ${matrix} -prefix ${aligned} ${part_mag}"

			3dAllineate -overwrite -base ${sbref} -final cubic \
				-1Dmatrix_apply ${matrix} \
				-prefix ${aligned} \
				${part_mag}
		fi
		part_mag=${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${ext}.nii.gz
		aligned=${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}_mcf_al${ext}.nii.gz

		echo "3dAllineate -overwrite -base ${sbref_ln} -final cubic -1Dmatrix_apply ${matrix} -prefix ${aligned} ${part_mag}"

		3dAllineate -overwrite -base ${sbref} -final cubic \
			-1Dmatrix_apply ${matrix} \
			-prefix ${aligned} \
			${part_mag} 

	done
done
