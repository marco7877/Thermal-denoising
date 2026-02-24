#!/bin/bash

module load afni/latest

source_directory=/bcbl/home/public/MarcoMotion/Rest_HighRes/S0

cd ${source_directory}

#list_subjects=("sub-001" "sub-002" "sub-003" "sub-004" "sub-005")
list_subjects=("sub-001")
list_methods=("vanilla" "mppca" "tmppca" "nordic" "hydra")
#list_tasks=("task-HABLA1200" "task-HABLA1700")
list_tasks=("task-REST")
#sub="${list_subjects[0]}"
#task="${list_tasks[0]}"
#method="${list_methods[0]}"
list_runs=( "run-1_")
#list_runs=( "")
target_directory=${source_directory}/imagesSPC
echo "Creating ${target_directory} if not exist"

mkdir -p ${target_directory}

export AFNI_NO_OBLIQUE_WARNING=YES

# opening AFNIb
echo "Starting loops"
#afni -niml -noplugins -yesplugouts &

for sub in "${list_subjects[@]}"; do
	echo "${sub}"
	for method in "${list_methods[@]}"; do
		echo "${method}"
		for task in "${list_tasks[@]}"; do
			echo "${task}"
			for run in "${list_runs[@]}"; do
				underlay=${sub}_ses-1_${task}_${run}echo-1_part-mag_masked_sbref
				cp /bcbl/home/public/MarcoMotion/Rest_HighRes/analysis_timeSeries/${underlay}.nii.gz ${source_directory}/${underlay}.nii.gz
				overlay=${sub}_ses-1_${task}_${run}S0meanSPC_part-mag_bold_${method}
				mask=${sub}_${method}_mask
				3dcalc -a ${underlay}.nii.gz -expr 'step(a)' -prefix ${mask}.nii.gz
				3dcalc -a ${sub}_ses-1_${task}_${run}S0mean_part-mag_bold_vanilla.nii.gz -b ${sub}_ses-1_${task}_${run}S0mean_part-mag_bold_${method}.nii.gz -expr '(b-a)/a' -prefix ${overlay}.nii.gz
				echo "Parameters:"
				echo "Underlay: ${underlay}"
				echo "Overlay: ${overlay}"

				sleep 2

				@chauffeur_afni \
					-pbar_posonly\
					-zerocolor Black\
					-ulay ${underlay}.nii.gz\
					-box_focus_slices ${overlay}.nii.gz\
					-olay ${overlay}.nii.gz\
					-olay_alpha Yes\
					-olay_boxed Yes\
					-no_cor \
					-opacity 9\
					-cbar yellow_to_cyan\
					-pbar_thr_alpha 0\
					-pbar_saveim ${target_directory}/${overlay}_pbar.png\
					-ulay_range 0% 120%\
					-func_range 0.05\
					-thr_olay 0\
					-label_mode 0\
					-prefix ${target_directory}/${overlay}\
					-save_ftype JPEG\
					-blowup 4\
					-montx 3\
					-monty 1\
					-montgap 0\
					-montcolor Black\
					-do_clean


				rm ${target_directory}/${overlay}_pbar*
				chmod 777 ${target_directory}/${overlay}*

				echo "Creating Montage"
				sleep 5
				mkdir -p  ${target_directory}/all
				2dcat\
					-gap -0\
					-gap_col '( 0 0 0 )'\
					-nx 1\
					-ny 2\
					-prefix ${target_directory}/all/${overlay}_all.jpg  ${target_directory}/${overlay}*.jpg
#					${target_directory}/${overlay}.axi.jpg ${target_directory}/${overlay}.sag.jpg 
				sleep 2
			done
		done
	done
done
#echo "Hit Enter to Quit AFNI"
#read pause
#plugout_drive  -com 'QUIT' \
	#	-quit   

