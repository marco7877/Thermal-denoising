#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=MPPCA_halves
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --mem=60G
#SBATCH --array=0-9
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt


module load MATLAB-IKUR


#list_task=( "task-HABLA1700" )
list_task=( "task-HABLA1200" "task-HABLA1700" )
#runs=( "1" "2" )
list_sub=( "001" "002" "003" "004" "005" )
method="mmpca"
split="True"
#list_sub=( "003" )
counter=0
for sub_label in "${list_sub[@]}" 
do
	for sub_task in "${list_task[@]}"
	do
				if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
					echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
					echo "Counter == ${counter}"
					break 2

				fi
				((counter++))
	done
done
sub_label=sub-$sub_label
echo "##################################################"
echo "############# SUBJECT :${sub_label}      #############"
echo "############# TASK :${sub_task}          #############"
echo "##################################################"
if [[ ${split} == "True" ]];then
	halves=( "1" "2" )
else
	halves=("")
fi
for half in "${halves[@]}"
do
mppca="/home/mflores/scripts/mppca.m"
temp_mppca=$(mktemp ~/mppcaXXXXXX.m)
origin=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}
output=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}

sed "s+SUBJECT+"${sub_label}"+;" $mppca>${temp_mppca}.e
sed "s+HALF+"${half}"+;" ${temp_mppca}.e>${temp_mppca}.ee
sed "s+TASK+"${sub_task}"+;" ${temp_mppca}.ee>$temp_mppca

matlab -batch "run('$temp_mppca');exit" 
rm ${temp_mppca}*
done

