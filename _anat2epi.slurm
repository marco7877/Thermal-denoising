#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=tsnr_unionMasks
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --array=1-30
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI/24.2.01-foss-2022a

echo "module load AFNI/24.3.06-foss-2022a"

echo "##################################################"
echo "$SLURM_JOB_NAME : $SLURM_ARRAY_JOB_ID, ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1700" )
list_subj=( "001" "002" "003" "004" "005" )
#runs=( 1 2 )
methods=( "nordic" "hydra" "vanilla" "tmmpca" "mppca" "nordic_magn" )
#methods=( "nordic" )
SE="False"
compute_mask=False
unionMask_tasks=True
skull_strip=False
allign_epi2anat=False
cp_files=True
counter=0
union_gm_wm=True
for task in "${tasks[@]}"
do
	for n_sub in "${list_subj[@]}" 
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
dir="/scratch/mflores/Resting_State"
preproc="func_preproc_${method}"
echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "############# TASK :${task}          #############"
echo "############# METHOD :${method}            #############"
echo "############# Single-echo : ${SE}            #############"
echo "##################################################"
subj=sub-${n_sub}

origin=${dir}/"${subj}/ses-1/func_preproc_vanilla"
origin_nat=${dir}/freesurfer/${subj}/SUMA
root=${dir}/"${subj}/ses-1/${preproc}"
matrix_suffix=a2e_${task}
analysis_folders=( "analysis" "analysis_timeSeries" )
echo analizing subject ${subj}

mask=${origin}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask.nii.gz
sbref=${origin}/${subj}_ses-1_${task}_echo-1_part-mag_sbref.nii.gz
t1=${origin_nat}/T1.nii.gz
skstr_t1=${origin_nat}/${subj}_ses-1_acq-uniclean-sk-str_T1w.nii.gz

masked_sbref=${root}/${subj}_ses-1_${task}_echo-1_part-mag_masked_sbref.nii.gz
alligned_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_alligned_mask.nii.gz
alligned_mask_union=${root}/${subj}_ses-1_${task}_echo-1_part-mag_alligned_mask_union.nii.gz
dilated_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_dilated_mask.nii.gz
non_gm_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm_mask.nii.gz
non_gm_mask_corrected=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm-corrected_mask.nii.gz
non_gm_mask_corrected_union=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm-corrected_mask_union.nii.gz
if [[ $SE == "True" ]];then
	target="${subj}_ses-1_${task}_SE_part-mag_bold_${method}"
else
	target="${subj}_ses-1_${task}_OC_part-mag_bold_${method}"
fi
#Computing GRAY MATTER and SUBCORTICAL masks from SUMA + FS

labeled_T1=${origin_nat}/aparc.a2009s+aseg.nii.gz

gm_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_gm-mask_T1w.nii.gz
slf_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_slf-mask_T1w.nii.gz
slb_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_slb-mask_T1w.nii.gz
sl_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_sl-mask_T1w.nii.gz
srf_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_srf-mask_T1w.nii.gz
srb_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_srb-mask_T1w.nii.gz
sr_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_sr-mask_T1w.nii.gz
sc_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_sc-mask_T1w.nii.gz
target_mask=${subj}_ses-1_${task}_${method}_echo-1_part-mag_alligned_mask.nii.gz
# gm
if [[ $compute_mask == True ]];then
	echo "making binary gm mask"
	3dcalc -a ${origin_nat}/lh.ribbon.nii.gz -b ${origin_nat}/rh.ribbon.nii.gz -expr '(a+b)' -prefix ${gm_mask} -overwrite

	# subcortical left hemisphere mask

	3dcalc -a ${labeled_T1} -expr 'within(a,9,13)' -prefix ${slf_mask} -overwrite

	3dcalc -a ${labeled_T1} -expr 'within(a,17,19)' -prefix ${slb_mask} -overwrite 

	3dcalc -a ${slf_mask} -b ${slb_mask} -expr '(a+b)' -prefix ${sl_mask} -overwrite

	rm ${slf_mask} ${slb_mask}

	# subcortical right hemisphere mask

	3dcalc -a ${labeled_T1} -expr 'within(a,48,51)' -prefix ${srf_mask} -overwrite

	3dcalc -a ${labeled_T1} -expr 'within(a,51,55)' -prefix ${srb_mask} -overwrite 

	3dcalc -a ${srf_mask} -b ${srb_mask} -expr '(a+b)' -prefix ${sr_mask} -overwrite

	rm ${srf_mask} ${srb_mask}

	# subcortical mask

	3dcalc -a ${sr_mask} -b ${sl_mask} -expr '(a+b)' -prefix ${sc_mask} -overwrite

	# target mask

	3dcalc -a ${gm_mask} -b ${sc_mask} -expr '(a+b)' -prefix ${origin_nat}/${subj}_ses-1_${task}_analysis-mask_T1w.nii.gz -overwrite


	3dmask_tool -dilate_inputs 1 -input ${origin_nat}/${subj}_ses-1_${task}_analysis-mask_T1w.nii.gz -prefix ${dilated_mask} -overwrite

	echo "masking ${sbref} with ${mask} as: ${masked_sbref}"

	3dcalc -a ${sbref} -b ${mask} -expr '(a*b)' -prefix ${masked_sbref} -overwrite

fi
if [[ $skull_strip == True ]];then
	echo "skull stripping ${t1} as: ${skstr_t1}"

	3dSkullStrip -input ${t1} -prefix ${skstr_t1} -overwrite
fi
if [[ $allign_epi2anat == True ]];then
	echo "Moving anat to epi"

	align_epi_anat.py -anat ${skstr_t1} -epi ${masked_sbref} -epi_base 0 -anat2epi -suffix ${matrix_suffix} -cost lpc -giant_move -anat_has_skull no -epi_strip None -volreg off -tshift off -output_dir ${root} -child_anat ${dilated_mask} -overwrite


	3dAllineate -base ${sbref} -1Dmatrix_apply ${root}/${subj}_ses-1_acq-uniclean-sk-str_T1w${matrix_suffix}_mat.aff12.1D -final NN -prefix ${alligned_mask} -input ${dilated_mask} -overwrite

	3dcalc -a ${mask} -b ${alligned_mask} -expr 'a*(a-b)' -prefix ${non_gm_mask} -overwrite

	3dmask_tool -input ${non_gm_mask} -dilate_inputs -1 1 -prefix ${non_gm_mask_corrected} -overwrite

	rm ${non_gm_mask}
fi
#rm ${gm_mask}
echo "##################################################"
echo " Using GM in the functional space to mask activation-"
echo "##################################################"

for directory in "${analysis_folders[@]}";do
	3dcalc -a ${dir}/${directory}/${target}.nii.gz -b ${alligned_mask} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_gm.nii.gz -overwrite
	3dcalc -a ${dir}/${directory}/${target}.nii.gz -b ${non_gm_mask_corrected} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_non-gm.nii.gz -overwrite
	cp ${alligned_mask} ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_gm-alligned_mask.nii.gz
	cp ${non_gm_mask_corrected} ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_wm-alligned_mask.nii.gz
	if [[ $unionMask_tasks == True ]];then
	3dcalc -a  ${dir}/${directory}/${target}_gm.nii.gz -b ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz  -expr '(a*b)' -prefix ${dir}/${directory}/${target}_gm_union.nii.gz -overwrite
	3dcalc -a  ${dir}/${directory}/${target}_non-gm.nii.gz -b ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz -expr '(a*b)' -prefix ${dir}/${directory}/${target}_non-gm_union.nii.gz -overwrite
	fi
	if [[ ${union_gm_wm} == True ]];then
		3dcalc -a ${alligned_mask} -b ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz  -expr '(a*b)' -prefix ${alligned_mask_union} -overwrite
	3dcalc -a ${non_gm_mask_corrected} -b ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz -expr '(a*b)' -prefix ${non_gm_mask_corrected_union} -overwrite

	fi
done


