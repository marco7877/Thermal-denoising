#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=MEICA_OC_1700_residuals
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --array=1-25
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load Anaconda3
source activate /data/mflores/TEDANA/



echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1700" )
list_subj=( "003" "001" "002" "004" "005" )
methods=(  "nordic" "tmmpca" "hydra" "mppca" "vanilla")
#method=( "nordic" )
#list_subj=( "01" "02" "03" "04" "06" "07" "08" "09" "11" "12" "14" "15" "16" "18" "19" "20" "21" )
#runs=( 1 2 )
#trees=( "tedana_orig" )
volume=""
counter=0
tedana="True"
t2smap="True"
residuals="True"
# Define echo times for each dataset 
echo_times="13.4 36.1 58.8 81.5"
#echo_times="11.2 28.1 45 61.9"
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

tree=meica


echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "############# TASK :${task}          #############"
echo "############# METHOD :${method}            #############"
echo "############# TREE :${tree}            #############"
echo "##################################################"

# define folders
subj=sub-${n_sub}
input=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}
output=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}/tedana_${task}_${tree}
nechoes=4
ext=${volume}.nii.gz

mkdir -p ${output}
mkdir -p ${output}/t2smap

if [[ ${residuals} = "True" ]];then
	part="residuals"
else
	part="bold"
fi

# Reference original volumes
part_mag=${input}/${subj}_ses-1_${task}_echo-1_part-mag_${part}_${method}_mcf_al${ext}
echo "My reference part_mag is: ${part_mag}"

# Define mask volume and create symbolic link

#mask_preproc=${input}/${subj}_ses-1_${task}_run-${run}_part-mag_desc-brain_mask.nii.gz
mask_preproc=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_vanilla/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz

echo "My reference mask is: ${mask_preproc}"

# make a list from 1 to 4 with 1 length digits
echo "List of echoes is ${list_echoes}"

# generate string with filenames of echoes for tedana
for n_echo in 2 3 4
do
	if [[ $residuals == "True" ]];then
		t2=" --t2smap ${output}bold${volume}/T2starmap.nii.gz "
	else
		t2=""
	fi
	mag=$part_mag" "${input}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_${part}_${method}_mcf_al${ext}
	part_mag=$mag
done

out_dir=${output}${part}${volume}
if [[ $tedana == "True" ]];then
	echo "##################################################"
	echo "tedana -d ${mag} -e ${echo_times} --out-dir ${out_dir} --tree ${tree} --verbose --debug True --mask ${mask_preproc}  --overwrite"
	echo "##################################################"

	tedana -d ${mag} -e ${echo_times} --out-dir ${out_dir} --tree ${tree} --verbose --debug --mask ${mask_preproc} --masktype none --fittype loglin ${t2}--overwrite

	rmdir --ignore-fail-on-non-empty ${out_dir}
fi
if [[ $t2smap == "True" ]];then
	echo "##################################################"
	echo "t2smap -d ${mag} -e ${echo_times} --out-dir ${out_dir}/t2smap${volume} --mask ${mask_preproc} --masktype none --fittype loglin --fitmode ts"
	echo "##################################################"
	t2smap -d ${mag} -e ${echo_times} --out-dir ${out_dir}/t2smap${volume} --mask ${mask_preproc} --masktype none --fittype loglin --fitmode ts 

	rmdir --ignore-fail-on-non-empty ${out_dir}/t2smap${volume}
fi


