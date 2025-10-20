#!/bin/bash 
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=Metrics_Rest_HighRes_all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --array=1-10
#SBATCH --error=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt

module load AFNI


echo "##################################################"
echo "$SLURM_JOB_NAME : $SLURM_ARRAY_JOB_ID, ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tsk="task-REST"
list_subj=( "001")
runs=( "1" "2" )
methods=("tmppca" "vanilla" "mppca" "nordic" "hydra")
#methods=( "vanilla" )
SE="False"
counter=0
TSNR="True"
T2="True"
R2="True"
halves="True"
residuals="False"
for n_sub in "${list_subj[@]}" 
do
	for run in "${runs[@]}"
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

task=${tsk}_run-${run}

echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "############# TASK :${task}_${run}          #############"
echo "############# METHOD :${method}            #############"
echo "##################################################"

if [[ $halves == "True" ]]; then
	volumes=( "" "1" "2" )
else
	volumes=( "" )
fi

for volume in "${volumes[@]}"; do

# define folders
subj=sub-${n_sub}
dir=/data/mflores/Rest_HighRes
origin=${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}
target="/scratch/mflores/Rest_HighRes/analysis"
target_T2="/scratch/mflores/Rest_HighRes/T2"
target_timeseries="/scratch/mflores/Rest_HighRes/analysis_timeSeries"
mkdir -p ${target} 
mkdir -p ${target_timeseries} 
mkdir -p ${target_T2}



part_mag=${origin}/desc-optcom_bold
target_mag=${target}/${subj}_ses-1_${task}_OC_part-mag_bold_${method}${volume}.nii.gz
part_magS=/scratch/mflores/Rest_HighRes/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}

if [[ $SE == "True" ]];then
	mkdir -p ${target_timeseries}
	origin=${dir}/${subj}/ses-1/func_preproc_${method}
	part_mag=${origin}/${subj}_ses-1_${task}_echo-2_part-mag_bold_${method}_mcf_al
	target_mag=${target}/${subj}_ses-1_${task}_SE_part-mag_bold_${method}.nii.gz
	cp ${part_mag}.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_SE_part-mag_bold_${method}.nii.gz
else
	cp ${part_mag}.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_OC_part-mag_bold_${method}${volume}.nii.gz
	if [[ $T2 == "True" ]];then
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap{volume}/T2starmap.nii.gz ${target_T2}/${subj}_ses-1_${task}_T2s_part-mag_bold_${method}${volume}.nii.gz
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/T2starmap.nii.gz ${target_T2}/${subj}_ses-1_${task}_T2_part-mag_bold_${method}${volume}.nii.gz
		echo "calculating signal percent change"

		3dTstat -mean -prefix ${target_T2}/${subj}_ses-1_${task}_T2mean_part-mag_bold_${method}${volume}.nii.gz ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap{volume}/T2starmap.nii.gz -overwrite 

		3dcalc -a ${target_T2}/${subj}_ses-1_${task}_T2mean_part-mag_bold_${method}${volume}.nii.gz -b ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap{volume}/T2starmap.nii.gz -expr '(a-b)/b' -prefix ${target_T2}/${subj}_ses-1_${task}_T2spc_part-mag_bold_${method}${volume}.nii.gz -overwrite


		#S0#
		cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/S0map.nii.gz ${target_T2}/${subj}_ses-1_${task}_S0_part-mag_bold_${method}${volume}.nii.gz
		echo "calculating signal percent change"
		3dcalc -a ${dir}/${subj}/ses-1/func_preproc_vanilla/tedana_${task}_meicabold${volume}/S0map.nii.gz -b ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/S0map.nii.gz -expr '(b-a)/a' -prefix ${target_T2}/${subj}_ses-1_${task}_S0spc_part-mag_bold_${method}${volume}.nii.gz -overwrite
	fi
	if [[ $R2 == "True" ]];then

		3dcalc -a ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicabold${volume}/t2smap{volume}/T2starmap.nii.gz -expr "(1/a)" -prefix ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz  -overwrite 

		3dTstat -mean -prefix ${target_T2}/${subj}_ses-1_${task}_R2mean_part-mag_bold_${method}${volume}.nii.gz ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz  -overwrite 

		3dcalc -a ${target_T2}/${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}.nii.gz -b ${target_T2}/${subj}_ses-1_${task}_R2mean_part-mag_bold_${method}${volume}.nii.gz -expr '(a-b)' -prefix ${target_T2}/${subj}_ses-1_${task}_R2spc_part-mag_bold_${method}${volume}.nii.gz -overwrite
	fi

fi

echo "My reference part_mag is: /${part_mag}"
mask=${dir}/${subj}/ses-1/func_preproc_vanilla/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz
ln_mask=${target}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz

if [[ $TSNR == "True" ]];then
	#OPTIMAL COMBINATION TSNR CALCULATION
	3dTstat -mean -prefix  ${part_magS}_mean.nii.gz ${part_mag}.nii.gz -overwrite


	3dTproject -polort 5 -prefix  ${part_magS}_dt.nii.gz -input ${part_mag}.nii.gz -overwrite


	3dcalc -a ${part_magS}_dt.nii.gz -b ${part_magS}_mean.nii.gz -expr 'a+b' -prefix ${part_magS}_dt.nii.gz -overwrite


	3dTstat -stdevNOD -prefix  ${part_magS}_std.nii.gz ${part_magS}_dt.nii.gz -overwrite


	3dcalc -a ${part_magS}_mean.nii.gz -b ${part_magS}_std.nii.gz -m ${mask} -expr 'm*(a/b)' -prefix ${part_magS}_tsnr.nii.gz -overwrite




	ln -s ${mask} ${ln_mask}

fi

if [[ $delta_t2 == "True" ]];then
	#OPTIMAL COMBINATION TSNR CALCULATION


	ln -s ${mask} ${ln_mask}
fi

if [[ $residuals == "True" ]];then

	echo "cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicaresiduals${volume}/desc-optcom_bold.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_residuals_part-mag_bold_${method}${volume}.nii.gz"
	cp ${dir}/${subj}/ses-1/func_preproc_${method}/tedana_${task}_meicaresiduals${volume}/desc-optcom_bold.nii.gz ${target_timeseries}/${subj}_ses-1_${task}_residuals_part-mag_bold_${method}${volume}.nii.gz
fi

done
