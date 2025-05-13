#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=trim_halves_nordic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:20:00
#SBATCH --mem=8G
#SBATCH --array=1-10
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI

echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1200" "task-HABLA1700" )
list_subj=( "003" "001" "002" "004" "005" )
#list_subj=( "01" "02" "03" "04" "06" "07" "08" "09" "11" "12" "14" "15" "16" "18" "19" "20" "21" )
#runs=( 1 2 )
methods=( "nordic" )
counter=0
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


#method=vanilla # vanilla, nordic, hydra 
preproc=func_preproc_${method} # better without the  
subj=sub-${n_sub}
# define folders

origin=/scratch/mflores/Resting_State/${subj}/ses-1/func
output=/scratch/mflores/Resting_State/${subj}/ses-1/${preproc}


trimvol=0 # number of volumes to trim due to steaty state magnetization

noise=3 # number of noise volumes
halves="False"
nechoes=4

# MAIN

mkdir -p ${output}


for parts in part-phase part-mag
do

	list_echoes=( "1" "2" "3" "4" )

	for n_echo in "${list_echoes[@]}"
	do
		echo "Trimming subject ${subj} nifti ${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz[${trimvol}..${vol}]"

		3dcalc -a "${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz[${trimvol}..$]" -expr 'a' \
			-prefix "${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz" -overwrite

		if [[ "${method}" == "vanilla" || "${method}" == "tmmpca" || "${method}" == "mppca" ]]
		then

			echo "3dinfo -nvi ${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz"

			volumes=$(3dinfo -nvi ${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz)

			vol=$(($volumes-$noise))

			echo "Original data ${volumes} volumes, but after removing noise volumes it has ${vol} volumes"
			3dcalc -a "${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz[${trimvol}..${vol}]" -expr 'a' \
				-prefix "${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz" -overwrite
			echo "3dcalc -a "${origin}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold.nii.gz[${trimvol}..${vol}]" -expr 'a' \
				-prefix "${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz" -overwrite"

			if [[ "${method}" == "vanilla" ]]
			then
				mv ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}.nii.gz
				echo "mv ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}.nii.gz"
			if [[ $halves == "True" ]]
			then
				epi_whole=${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}.nii.gz
				n_epi_whole=$(3dinfo -nvi ${epi_whole})
				epi_half=$(($n_epi_whole/2))
				echo "Triming ${epi_whole} in halves"
				3dcalc -a "${epi_whole}[0..${epi_half}]" -expr "a" -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}1.nii.gz -overwrite

				3dcalc -a "${epi_whole}[${epi_half}..$]" -expr "a" -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}2.nii.gz -overwrite

			fi
			fi
			if [[ $halves == "True" ]]
			then
				epi_whole=${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd.nii.gz
				n_epi_whole=$(3dinfo -nvi ${epi_whole})
				epi_half=$(($n_epi_whole/2))
				echo "Triming ${epi_whole} in halves"
				3dcalc -a "${epi_whole}[0..${epi_half}]" -expr "a" -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd1.nii.gz -overwrite

				3dcalc -a "${epi_whole}[${epi_half}..$]" -expr "a" -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_${parts}_bold_${method}_dsd2.nii.gz -overwrite

			fi
		fi

	done
done


