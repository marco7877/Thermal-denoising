#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=True_Resting_State_OC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --array=1-50
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI


echo "##################################################"
echo "$SLURM_JOB_NAME : $SLURM_ARRAY_JOB_ID, ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1200" "task-HABLA1700" )
list_subj=( "001" "002" "003" "004" "005" )
#runs=( 1 2 )
methods=( "vanilla" "hydra" "nordic" "tmmpca" "mppca" )
#methods=( "nordic" )
SE="False"
counter=0
TSNR="True"
T2="False"
R2="False"
volume=""
for n_sub in "${list_subj[@]}" 
do
	for task in "${tasks[@]}"
	do
		for method in "${methods[@]}"
		do
			((counter++))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 3

			fi
		done
	done
done



echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "############# TASK :${task}          #############"
echo "############# METHOD :${method}            #############"
echo "##################################################"

# define folders
subj=sub-${n_sub}
dir=/scratch/mflores/Resting_State
origin=${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}
target="${dir}/analysis"
target_T2=${dir}/T2
target_timeseries="${dir}/analysis_timeSeries"
mkdir -p ${target} 
mkdir -p ${target_timeseries} 
mkdir -p ${target_T2}



part_mag=${origin}/desc-optcom_bold
target_mag=${target}/${subj}_ses-1_${task}_OC_part-mag_bold_${method}${volume}.nii.gz

if [[ $SE == "True" ]];then
	mkdir -p ${target_timeseries}
	origin=${dir}/${subj}/ses-1/func_preproc_${method}
	part_mag=${origin}/${subj}_ses-1_${task}_echo-2_part-mag_bold_${method}_mcf_al
	target_mag=${target}/${subj}_ses-1_${task}_SE_part-mag_bold_${method}.nii.gz
	cp ${part_mag}.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_SE_part-mag_bold_${method}.nii.gz
else
	cp ${part_mag}.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_OC_part-mag_bold_${method}${volume}.nii.gz
	if [[ $T2 == "True" ]];then
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap/T2starmap.nii.gz ${target_T2}/${subj}_ses-1_${task}_T2s_part-mag_bold_${method}${volume}.nii.gz
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/T2starmap.nii.gz ${target_T2}/${subj}_ses-1_${task}_T2_part-mag_bold_${method}${volume}.nii.gz
		echo "calculating signal percent change"

		3dTstat -mean -prefix ${target_T2}/${subj}_ses-1_${task}_T2mean_part-mag_bold_${method}${volume}.nii.gz ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap/T2starmap.nii.gz -overwrite 

		3dcalc -a ${target_T2}/${subj}_ses-1_${task}_T2mean_part-mag_bold_${method}${volume}.nii.gz -b ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap/T2starmap.nii.gz -expr '(a-b)/b' -prefix ${target_T2}/${subj}_ses-1_${task}_T2spc_part-mag_bold_${method}${volume}.nii.gz -overwrite


		#S0#
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/S0map.nii.gz ${target_T2}/${subj}_ses-1_${task}_S0_part-mag_bold_${method}${volume}.nii.gz
		echo "calculating signal percent change"
		3dcalc -a ${dir}/${subj}/ses-1/func_preproc_vanilla/tedana_${task}_meicabold${volume}/S0map.nii.gz -b ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/S0map.nii.gz -expr '(b-a)/a' -prefix ${target_T2}/${subj}_ses-1_${task}_S0spc_part-mag_bold_${method}${volume}.nii.gz -overwrite
		if [[ $R2 == "True" ]];then

			3dcalc -a ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap/T2starmap.nii.gz -expr "(1/a)" -prefix ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz  -overwrite 

			3dTstat -mean -prefix ${target_T2}/${subj}_ses-1_${task}_R2mean_part-mag_bold_${method}${volume}.nii.gz ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz  -overwrite 

			3dcalc -a ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz -b ${target_T2}/${subj}_ses-1_${task}_R2mean_part-mag_bold_${method}${volume}.nii.gz -expr '(a-b)' -prefix ${target_T2}/${subj}_ses-1_${task}_R2spc_part-mag_bold_${method}${volume}.nii.gz -overwrite
		fi

	fi

fi
echo "My reference part_mag is: /${part_mag}"
mask=${dir}/${subj}/ses-1/func_preproc_vanilla/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz
ln_mask=${target}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz

if [[ $TSNR == "True" ]];then
	#OPTIMAL COMBINATION TSNR CALCULATION
	3dTstat -mean -prefix  ${part_mag}_mean.nii.gz ${part_mag}.nii.gz -overwrite


	3dTproject -polort 5 -prefix  ${part_mag}_dt.nii.gz -input ${part_mag}.nii.gz -overwrite


	3dcalc -a ${part_mag}_dt.nii.gz -b ${part_mag}_mean.nii.gz -expr 'a+b' -prefix ${part_mag}_dt.nii.gz -overwrite


	3dTstat -stdevNOD -prefix  ${part_mag}_std.nii.gz ${part_mag}_dt.nii.gz -overwrite


	3dcalc -a ${part_mag}_mean.nii.gz -b ${part_mag}_std.nii.gz -m ${mask} -expr 'm*(a/b)' -prefix ${part_mag}_tsnr.nii.gz -overwrite


	cp ${part_mag}_tsnr.nii.gz ${target_mag}


	ln -s ${mask} ${ln_mask}

fi

if [[ $delta_t2 == "True" ]];then
	#OPTIMAL COMBINATION TSNR CALCULATION
	cp ${part_mag}_tsnr.nii.gz ${target_mag}


	ln -s ${mask} ${ln_mask}

	if
