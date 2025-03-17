#!/bin/bash
#$ -m be
#$ -m be
#$ -M m.flores@bcbl.eu
#$ -S /bin/bash



#list_subj=$( count -digits 1 1 5 )
tasks=( task-WORD task-SENT task-SYLAB )
runs=( 1 )

list_subj1=( 06 08 09 11 12 13 14 15 16 17 18 19 20 21 )
list_subj=${list_subj1[0]}
echo "List of subjects is ${list_subj}"

for n_sub in ${list_subj[@]}
do
	for task in ${tasks[@]}
	do
		for run in ${runs[0]}
		do
			error_txt=/scratch/mflores/HABLA/logs/sub${n_sub}_${task}_${run}_ME-ICA.txt

			if [[ -e ${error_txt} ]]; then

				rm ${error_txt}

			fi

#			echo "sbatch  --job-name=sub${n_sub}_${task}_${run}_ME-ICA --output=${error_txt} --error=${error_txt} /scratch/mflores/scripts/_ME-ICA.slurm ${n_sub} ${task} ${run}" 
			qsub -q long.q -N sub${n_sub}_${task}_${run}_OC -o ${error_txt} -e ${error_txt} /bcbl/home/public/MarcoMotion/scripts/HABLA_SPiN/MEICA.sh ${n_sub} ${task} ${run} 
		done
	done
done
