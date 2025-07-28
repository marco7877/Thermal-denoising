#!/bin/bash -xv
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=HYDRA_Rest_HighRes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mem=90G
#SBATCH --array=1-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Rest_HighRes/logs/%x_%a_%j.txt



module load AFNI/24.3.06-foss-2022a
module load MATLAB-IKUR

echo "##################################################"
echo "HYDRA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"

list_task=( "task-REST" )
runs=( "1" "2" )
list_subj=( "001")
counter=0
make_residuals="True"
halves="True"
for n_sub in "${list_subj[@]}" 
do
	for tsk in "${list_task[@]}"
	do
		for run in "${runs[@]}"
		do
		((counter++))
		if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
			echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
			break 3

		fi
	done
done
done
method=hydra
preproc=func_preproc_${method}


subj=sub-${n_sub}
task=${tsk}_run-${run}
#task=${tsk}
repo=/home/mflores/scripts
origin=/data/mflores/Rest_HighRes/${subj}/ses-1/func_preproc_vanilla
output=/scratch/mflores/Rest_HighRes/${subj}/ses-1/${preproc}
noise=3 # number of noise volumes
save_gfactor=1 # save gfactor map = 1 (saves gfactor and runs nordic); = 2 (saves gfactor and stops)  
nechoes=3
# MAIN

#Check if output durectory exists

echo "##################################################"
echo "############# SUBJECT :${subj}      #############"
echo "############# TASK :${task}          #############"
echo "##################################################"


mkdir -p ${output}/tmp


echo "*************************************************************"
# Reference outputal volumnes
part_mag=/data/mflores/Rest_HighRes/${subj}/ses-1/${preproc}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_dsd.nii.gz

echo "My reference part_mag is: ${part_mag}"

part_phase=/data/mflores/Rest_HighRes/${subj}/ses-1/${preproc}/${subj}_ses-1_${task}_echo-1_part-phase_bold_${method}_dsd.nii.gz

echo "My reference part_phase is: ${part_phase}"



if [[ $halves == "True" ]]; then
	extentions=( "" "1" "2" )
else
	extentions=( "" )
fi

# make a list from 1 to 4 with 1 length digits

list_echoes=($(seq 1 ${nechoes}))

echo "List of echoes is ${list_echoes}"
list_aggregate=($(seq 2 ${nechoes}))

for n_echo in "${list_aggregate[@]}"

do

	mag=$part_mag" "/data/mflores/Rest_HighRes/${subj}/ses-1/${preproc}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}_dsd.nii.gz

	phase=$part_phase" "/data/mflores/Rest_HighRes/${subj}/ses-1/${preproc}/${subj}_ses-1_${task}_echo-${n_echo}_part-phase_bold_${method}_dsd.nii.gz

	part_mag=$mag

	part_phase=$phase

done


echo "Magnitude volumes to concatenate are: $part_mag"

echo "Phase volumes to concatenate are: $part_phase"
3dZcat $part_mag -prefix ${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd.nii.gz -overwrite

3dZcat $part_phase -prefix ${output}/${subj}_ses-1_${task}_echoes_part-phase_bold_${method}_dsd.nii.gz -overwrite


if [[ $halves == "True" ]]; then

	volumes=$(3dinfo -nvi ${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd.nii.gz)
	noise_vol=$(($volumes-$(($noise-1))))
	signal_vol=$(($noise_vol-1))
	vol_half=$(($noise_vol/2))
	vol_half_end=$(($vol_half-1))
	echo "Unsplitted time series has ${volumes} volumes"
	echo "After taking out noise scans we have ${noise_vol} volumes"
	echo "Half of the data set is ${vol_half}"

		part_mag="${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd"

		part_mag_half1="${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd1"

		part_mag_half2="${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd2"

		part_mag_out1="${subj}_ses-1_${task}_echoes_part-mag_bold_${method}N1"

		part_mag_out2="${subj}_ses-1_${task}_echoes_part-mag_bold_${method}N2"

		part_phase="${output}/${subj}_ses-1_${task}_echoes_part-phase_bold_${method}_dsd"

		part_phase_half1="${output}/${subj}_ses-1_${task}_echoes_part-phase_bold_${method}_dsd1"

		part_phase_half2="${output}/${subj}_ses-1_${task}_echoes_part-phase_bold_${method}_dsd2"

		noise_mag="${output}/${subj}_ses-1_${task}_echoes_part-mag_noise_${method}"

		noise_phase="${output}/${subj}_ses-1_${task}_echoes_part-phase_noise_${method}"

		echo "Cutting volumes in halves"

		echo "Unsplitted time series has ${volumes} volumes"
		echo "After taking out noise scans we have ${noise_vol} volumes"
		echo "Half of the data set is ${vol_half}"

		echo  "3dcalc  -a "${part_mag}.nii.gz[${noise_vol}..$]" -expr "a" -prefix ${noise_mag}.nii.gz -overwrite"
		3dcalc -a "${part_mag}.nii.gz[${noise_vol}..$]" -expr "a" -prefix ${noise_mag}.nii.gz -overwrite
		echo "3dcalc -a "${part_phase}.nii.gz[${noise_vol}..$]" -expr "a" -prefix ${noise_phase}.nii.gz -overwrite"
		3dcalc -a "${part_phase}.nii.gz[${noise_vol}..$]" -expr "a" -prefix ${noise_phase}.nii.gz -overwrite


		3dcalc -a "${part_mag}.nii.gz[0..${vol_half_end}]" -expr "a" -prefix ${part_mag_half1}u.nii.gz -overwrite

		3dTcat -prefix ${part_mag_half1}.nii.gz ${part_mag_half1}u.nii.gz "${noise_mag}.nii.gz[1]" -overwrite

		rm ${part_mag_half1}u.nii.gz
		
		3dcalc -a "${part_mag}.nii.gz[${vol_half}..${signal_vol}]" -expr "a" -prefix ${part_mag_half2}u.nii.gz -overwrite

		3dTcat -prefix ${part_mag_half2}.nii.gz ${part_mag_half2}u.nii.gz "${noise_mag}.nii.gz[2]" -overwrite
		
		rm ${part_mag_half2}u.nii.gz
		
		3dcalc -a "${part_phase}.nii.gz[0..${vol_half_end}]" -expr "a" -prefix ${part_phase_half1}u.nii.gz -overwrite

		3dTcat -prefix ${part_phase_half1}.nii.gz ${part_phase_half1}u.nii.gz "${noise_phase}.nii.gz[1]" -overwrite
		
		rm ${part_phase_half1}u.nii.gz

		3dcalc -a "${part_phase}.nii.gz[${vol_half}..${signal_vol}]" -expr "a" -prefix ${part_phase_half2}u.nii.gz -overwrite

		3dTcat -prefix ${part_phase_half2}.nii.gz ${part_phase_half2}u.nii.gz "${noise_phase}.nii.gz[2]" -overwrite
		
		rm ${part_phase_half2}u.nii.gz
		
		mags=(${part_mag_half1} ${part_mag_half2})
		
		phases=(${part_phase_half1} ${part_phase_half2})
		
		outs=(${part_mag_out1} ${part_mag_out2})
fi
for ext in "${extentions[@]}"; do
if [[ $ext == "" ]]; then
	noise=3
else
	noise=1
fi
	temp_nordic=$(mktemp ${output}/tmp/nordicXXXXXX.m)
	
	echo "*************************************************************"
	echo "addpath(genpath('/home/mflores/nordic'));cd '${repo}';ARG.temporal_phase=1;ARG.noise_volume_last=${noise};ARG.phase_filter_width=10;ARG.save_add_info=1;ARG.save_gfactor_map=${save_gfactor};ARG.DIROUT='${output}/';NIFTI_NORDIC('${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd${ext}.nii.gz','${output}/${subj}_ses-1_${task}_echoes_part-phase_bold_${method}_dsd${ext}.nii.gz','${subj}_ses-1_${task}_echoes_part-mag_bold_${method}${ext}',ARG);" > $temp_nordic

	matlab -batch "run('$temp_nordic');exit" 


	gzip -f ${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}${ext}.nii


	z_orig=$(3dinfo -nk ${origin}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}_dsd${ext}.nii.gz
	)
	start=0

	finish=$(($z_orig-1))

	for n_echoe in "${list_echoes[@]}"

	do

		echo "*************************************************************"
		echo "3dZcutup -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz -keep ${start} ${finish} ${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}_dsd${ext}.nii.gz -overwrite"

		3dZcutup -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz -keep ${start} ${finish} ${output}/${subj}_ses-1_${task}_echoes_part-mag_bold_${method}${ext}.nii.gz -overwrite


		start=$(($start+$z_orig))

		finish=$(($finish+$z_orig))
	done
	
	ATR=$( 3dAttribute IJK_TO_DICOM_REAL ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}${ext}.nii.gz )

	volumes=$(3dinfo -nvi ${output}/${subj}_ses-1_${task}_echo-1_part-mag_bold_${method}${ext}.nii.gz)

	vol=$(($volumes-$noise))

	echo "I have original $volumes  volumes but after removing noise we have: $vol"

	for n_echoe in "${list_echoes[@]}"

	do

		3drefit -atrfloat IJK_TO_DICOM_REAL "${ATR}" ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz

		# Removing last noise volumes

		echo "removing noise to ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz"

		3dcalc -a "${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz[0..${vol}]" -expr 'a' -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz -overwrite
	if [[ $make_residualds == "True" ]]; then
		echo "Making residuals between vanilla and denoised volumes"
		3dcalc -a ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_${method}${ext}.nii.gz -b ${origin}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_bold_vanilla${ext}.nii.gz -expr '(b-a)' -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echoe}_part-mag_residuals_${method}${ext}.nii.gz
	fi
done
done


#rm -rf ${output}tmp
