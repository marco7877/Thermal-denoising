#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=OC_weights
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=30G
#SBATCH --array=1-10
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI


echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1700" "task-HABLA1200" )
list_subj=( "003" "001" "002" "004" "005" )
method=("vanilla")
residuals="True"
for n_sub in "${list_subj[@]}" 
do
	for task in "${tasks[@]}"
	do
			((counter++))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 2

			fi
	done
done


echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "############# TASK :${task}          #############"
echo "############# METHOD :${method}            #############"
echo "############# RESIDUALS :${residuals}            #############"
echo "##################################################"

# define folders
subj=sub-${n_sub}
input=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_vanilla/tedana_${task}_meicabold

echoes=4
if [[ ! -e ${input}/T2starmap.nii.gz  ]]; then
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "Niifti ${input}/T2starmap.nii.gz does not exist"
	echo "Try runing tedana before runing this script"
	echo "Exiting ..."
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	exit
fi

echo "-------------------------------------------------------------"
echo "Creating echo weights!"
echo "-------------------------------------------------------------"

if [[ ${residuals} = "True" ]];then
	part="residuals"
else
	part="bold"
fi

if [[ ${task} = "task-HABLA1700" ]]
then
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(13.4*(exp(13.4/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-1_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(13.4*(exp(13.4/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-1_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(36.1*(exp(36.1/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-2_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(36.1*(exp(36.1/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-2_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(58.8*(exp(58.8/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-3_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(58.8*(exp(58.8/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-3_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(81.5*(exp(81.5/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-4_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(81.5*(exp(81.5/a)*-1))/((13.4*(exp(13.4/a)*-1)+(36.1*(exp(36.1/a)*-1)+(58.8*(exp(58.8/a)*-1)+(81.5*(exp(81.5/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-4_weight.nii.gz -overwrite

elif [[ ${task} = "task-HABLA1200" ]]
then
	echo_times="11.2 28.1 45 61.9"
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(11.2*(exp(11.2/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-1_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(11.2*(exp(11.2/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-1_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(28.1*(exp(28.1/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-2_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(28.1*(exp(28.1/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-2_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(45*(exp(45/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-3_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(45*(exp(45/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-3_weight.nii.gz -overwrite
	
	echo "-------------------------------------------------------------"
	echo "3dcalc -a ${input}/T2starmap.nii.gz -expr '(61.9*(exp(61.9/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-4_weight.nii.gz -overwrite"
	3dcalc -a ${input}/T2starmap.nii.gz -expr '(61.9*(exp(61.9/a)*-1))/((11.2*(exp(11.2/a)*-1)+(28.1*(exp(28.1/a)*-1)+(45*(exp(45/a)*-1)+(61.9*(exp(61.9/a)*-1))))))' -prefix ${input}/${subj}_${task}_OC_echo-4_weight.nii.gz -overwrite
fi
