#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=MEICA_OC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --array=1-40
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load Anaconda3
source activate /data/mflores/TEDANA/

#module load AFNI


echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1700" "task-HABLA1200" )
list_subj=( "003" "001" "002" "004" "005" )
methods=(  "nordic" "tmmpca" "hydra" "mppca" )

volume=""
counter=0
tedana="True"
t2smap="True"
residuals="False"
# ir residuals == true tedana and t2smap will not be useful
# Define echo times for each dataset 

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
if [[ ${task} = "task-HABLA1700" ]];then
	echo_times="13.4 36.1 58.8 81.5"
elif [[ ${task} = "task-HABLA1200" ]];then
	echo_times="11.2 28.1 45 61.9"
fi


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
	mag=$part_mag" "${input}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_${part}_${method}_mcf_al${ext}
	part_mag=$mag
done

out_dir=${output}${part}${volume}
if [[ $residuals == "True" ]];then
input_weights=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_vanilla/tedana_${task}_${tree}bold${volume}
	echo "##################################################"
	echo "Making OC from residuals"
	echo "##################################################"
part_mag1=${input}/${subj}_ses-1_${task}_echo-1_part-mag_${part}_${method}_mcf_al${ext}
part_mag2=${input}/${subj}_ses-1_${task}_echo-2_part-mag_${part}_${method}_mcf_al${ext}
part_mag3=${input}/${subj}_ses-1_${task}_echo-3_part-mag_${part}_${method}_mcf_al${ext}
part_mag4=${input}/${subj}_ses-1_${task}_echo-4_part-mag_${part}_${method}_mcf_al${ext}
weight_mag1=${input_weights}/${subj}_${task}_OC_echo-1_weight.nii.gz
weight_mag2=${input_weights}/${subj}_${task}_OC_echo-2_weight.nii.gz
weight_mag3=${input_weights}/${subj}_${task}_OC_echo-3_weight.nii.gz
weight_mag4=${input_weights}/${subj}_${task}_OC_echo-4_weight.nii.gz

3dcalc -a ${part_mag1} -b ${part_mag2} -c ${part_mag3} -d ${part_mag4} -e ${weight_mag1} -f ${weight_mag2} -g ${weight_mag3} -h ${weight_mag4} -i ${mask_preproc}  -expr '(a*e+b*f+c*g+d*h)*i' -prefix /scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}/tedana_${task}_${tree}${part}${volume}/desc-optcom_bold.nii.gz -overwrite 
	else

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
fi
