#!/bin/bash -xv
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=MPPCA_MEVOTLOC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH --mem=60G
#SBATCH --array=0-7
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/MEVOTLOC/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/MEVOTLOC/logs/%x_%a_%j.txt


module load MATLAB-IKUR
module load AFNI

PS4='${LINENO}: '
#list_task=( "task-HABLA1700" )
list_task=( "task-fLoc" )
runs=( "1" "2" "3" "4" "5" "6" "7" )
list_sub=( "001" )
method="mppca"
split="True"
make_residuals="True"
nechoes=3
#list_sub=( "003" )
counter=0
for sub_label in "${list_sub[@]}" 
do
	for tsk in "${list_task[@]}"
	do
		for run in "${runs[@]}"
		do
				if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
					echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
					echo "Counter == ${counter}"
					break 3

				fi
				((counter++))
	done
done
done
subj=sub-$sub_label

task=${tsk}_run-${run}
echo "##################################################"
echo "############# SUBJECT :${subj}      #############"
echo "############# TASK :${task}          #############"
echo "##################################################"
if [[ ${split} == "True" ]];then
	halves=( "1" "2" "" )
else
	halves=("")
fi
for half in "${halves[@]}"
do
directory=/data/mflores/MEVOTLOC
mppca="/home/mflores/Thermal-denoising/mppca.m"
temp_mppca=$(mktemp ~/mppcaXXXXXX.m)
origin=${directory}/${subj}/ses-1/func_preproc_vanilla
output=/scratch/mflores/MEVOTLOC/${subj}/ses-1/func_preproc_${method}

sed "s+SUBJECT+"${subj}"+;" $mppca>${temp_mppca}.e
sed "s+DIRECTORY+"${directory}"+;" ${temp_mppca}.e>${temp_mppca}.i
sed "s+HALF+"${half}"+;" ${temp_mppca}.i>${temp_mppca}.ee
sed "s+TASK+"${task}"+;" ${temp_mppca}.ee>$temp_mppca

matlab -batch "run('$temp_mppca');exit" 
rm ${temp_mppca}*
	if [[ $make_residuals == "True" ]]; then
		echoes=($(seq 1 ${nechoes}))
		echo "Making residuals between vanilla and denoised volumes"
		for n_echo in "${echoes[@]}"; do

			3dcalc -a ${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_${method}${half}.nii.gz -b ${origin}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_bold_vanilla${half}.nii.gz -expr '(b-a)' -prefix ${output}/${subj}_ses-1_${task}_echo-${n_echo}_part-mag_residuals_${method}${half}.nii.gz -overwrite
		done
	fi

done

