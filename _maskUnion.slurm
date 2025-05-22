#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=goodSignal_mask_union
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --mem=30G
#SBATCH --array=1-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI/24.3.06-foss-2022a

${SLURM_ARRAY_TASK_ID}
n_subs=( "sub-001" "sub-002" "sub-003" "sub-004" "sub-005" )
methods=( "vanilla" "nordic" "mppca" "nordic_magn" "tmmpca" )
tasks=( "task-HABLA1200" "task-HABLA1700" )
counter=0
volume=""
for sub in "${n_subs[@]}" 
do
	for task in "${tasks[@]}"
	do
			(( counter++ ))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 2

	fi
	# echo " subject: ${sub} , task: ${task}, counter: ${counter}"
	done
done



echo"########################################"
echo" subject ${sub}"
echo" Creating union mask for ${task}"
echo"########################################"


dir=/scratch/mflores/Resting_State/${sub}/ses-1/func_preproc

masks=/scratch/mflores/Resting_State/analysis/${sub}_ses-1_${task}_echo-1_part-mag_brain_mask-bucket.nii.gz
mask_union=/scratch/mflores/Resting_State/analysis/${sub}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz
mask_union_ts=/scratch/mflores/Resting_State/analysis_timeSeries/${sub}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz
for method in "${methods[@]}";do
	method_mask=${dir}_${method}/tedana_${task}_meicabold${volume}/desc-adaptiveGoodSignal_mask.nii.gz
	method_mask_binary=${dir}_${method}/tedana_${task}_meicabold${volume}/desc-adaptiveGoodSignal_mask_binary.nii.gz
	3dcalc -a ${method_mask} -expr 'astep(a,1)' -overwrite -prefix ${method_mask_binary}
	3dbucket -aglueto ${masks} -fbuc ${method_mask_binary}
##	masks=${masks}" "${method_mask_binary}
done
3dmask_tool -input ${masks} -union -prefix ${mask_union} -overwrite
 cp ${mask_union} ${mask_union_ts} 

