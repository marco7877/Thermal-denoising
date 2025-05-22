#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=epi2anat_union-gm_half2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=6G
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --array=1-10
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI/24.2.01-foss-2022a

echo "module load AFNI/24.3.06-foss-2022a"

echo "##################################################"
echo "$SLURM_JOB_NAME : $SLURM_ARRAY_JOB_ID, ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
tasks=( "task-HABLA1700" "task-HABLA1200" )
list_subj=( "001" "002" "003" "004" "005" )
#runs=( 1 2 )
#methods=( "nordic" "hydra" "vanilla" "tmmpca" "mppca" )
methods=( "nordic" )
SE="False"
volume="2"
compute_mask="False"
skull_strip="False"
allign_epi2anat="False"
counter=0
union_gm_wm="True"
residuals="False"
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
sbref=/data/HABLA/Resting_State/${subj}/ses-1/func/${subj}_ses-1_${task}_echo-1_part-mag_sbref.nii.gz
t1=${origin_nat}/T1.nii.gz
skstr_t1=${origin_nat}/${subj}_ses-1_acq-uniclean-sk-str_T1w.nii.gz

masked_sbref=/scratch/mflores/Resting_State/${subj}/ses-1/func_preproc_${method}/${subj}_ses-1_${task}_echo-1_part-mag_masked_sbref.nii.gz
alligned_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_alligned_mask.nii.gz
alligned_mask_union=${root}/${subj}_ses-1_${task}_echo-1_part-mag_alligned_mask_union.nii.gz
dilated_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_dilated_mask.nii.gz
non_gm_mask=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm_mask.nii.gz
non_gm_mask_corrected=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm-corrected_mask.nii.gz
non_gm_mask_corrected_union=${root}/${subj}_ses-1_${task}_echo-1_part-mag_non-gm-corrected_mask_union.nii.gz
if [[ $SE == "True" ]];then
	target="${subj}_ses-1_${task}_SE_part-mag_bold_${method}"
else
	target="${subj}_ses-1_${task}_OC_part-mag_bold_${method}${volume}"
	target_T2="${subj}_ses-1_${task}_T2_part-mag_bold_${method}${volume}"
	target_R2="${subj}_ses-1_${task}_R2_part-mag_bold_${method}${volume}"
	target_deltaR2="${subj}_ses-1_${task}_R2spc_part-mag_bold_${method}${volume}"
	target_S0="${subj}_ses-1_${task}_S0_part-mag_bold_${method}${volume}"
	target_residuals="${subj}_ses-1_${task}_residuals_part-mag_bold_${method}${volume}"

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
left_cerebellum=${origin_nat}/${subj}_ses-1_${task}_${method}_lcerebellum-mask_T1w.nii.gz
right_cerebellum=${origin_nat}/${subj}_ses-1_${task}_${method}_rcerebellum-mask_T1w.nii.gz
brainsteam=${origin_nat}/${subj}_ses-1_${task}_${method}_brainsteam-mask_T1w.nii.gz
CSF=${origin_nat}/${subj}_ses-1_${task}_${method}_csf-mask_T1w.nii.gz
areas_not_interest=${origin_nat}/${subj}_ses-1_${task}_${method}_notbrain-mask_T1w.nii.gz
areas_not_interest_epi=${origin_nat}/${subj}_ses-1_${task}_${method}_notbrain-mask_epi.nii.gz
wm_mask_seg=${origin_nat}/wm.nii.gz
wm_mask_seg_binary=${origin_nat}/wm_binary.nii.gz
wm_mask_tmp=${origin_nat}/${subj}_ses-1_${task}_${method}_wm-mask_tmpT1w.nii.gz
wm_mask=${origin_nat}/${subj}_ses-1_${task}_${method}_wm-mask_T1w.nii.gz


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

	# cerebellum, CSF and brainsteam mask

	3dcalc -a ${labeled_T1} -expr 'within(a,46,47)' -prefix ${left_cerebellum} -overwrite

	3dcalc -a ${labeled_T1} -expr 'within(a,7,8)' -prefix ${right_cerebellum} -overwrite

	3dcalc -a ${labeled_T1} -expr 'equals(a,16)' -prefix ${brainsteam} -overwrite

	3dcalc -a ${labeled_T1} -expr 'equals(a,24)' -prefix ${CSF} -overwrite

	3dcalc -a ${left_cerebellum} -b ${right_cerebellum} -c ${CSF} -d ${brainsteam} -expr '(a+b+c+d)' -prefix ${areas_not_interest} -overwrite

	3dcalc -a ${wm_mask_seg} -expr 'ispositive(a)' -prefix ${wm_mask_seg_binary} -overwrite

	3dcalc -a ${wm_mask_seg_binary} -b ${areas_not_interest} -expr '(a-b)' -prefix ${wm_mask_tmp} -overwrite
	3dcalc -a ${wm_mask_tmp} -expr 'ispositive(a)' -prefix ${wm_mask} -overwrite
	

	rm ${left_cerebellum} ${right_cerebellum} ${brainsteam} ${CSF}


	# target mask

	3dcalc -a ${gm_mask} -b ${sc_mask} -expr '(a+b)' -prefix ${origin_nat}/${subj}_ses-1_${task}_${method}_analysis-mask_T1w.nii.gz -overwrite


	3dmask_tool -input ${origin_nat}/${subj}_ses-1_${task}_${method}_analysis-mask_T1w.nii.gz -prefix ${dilated_mask} -overwrite

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

	3dAllineate -base ${sbref} -1Dmatrix_apply ${root}/${subj}_ses-1_acq-uniclean-sk-str_T1w${matrix_suffix}_mat.aff12.1D -final NN -prefix ${areas_not_interest_epi} -input ${areas_not_interest} -overwrite

	3dAllineate -base ${sbref} -1Dmatrix_apply ${root}/${subj}_ses-1_acq-uniclean-sk-str_T1w${matrix_suffix}_mat.aff12.1D -final NN -prefix ${non_gm_mask_corrected} -input ${wm_mask} -overwrite
	

#	3dmask_tool -input ${non_gm_mask} -prefix ${non_gm_mask_corrected} -overwrite

#	rm ${non_gm_mask}
fi
#rm ${gm_mask}
echo "##################################################"
echo " Using GM in the functional space to mask tSNR maps-"
echo "##################################################"

for directory in "${analysis_folders[@]}";do
	3dcalc -a ${dir}/${directory}/${target}.nii.gz -b ${alligned_mask} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_gm.nii.gz -overwrite
	3dcalc -a ${dir}/${directory}/${target}.nii.gz -b ${non_gm_mask_corrected} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_non-gm.nii.gz -overwrite
	cp ${alligned_mask} ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_gm-alligned_mask.nii.gz
	cp ${non_gm_mask_corrected} ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_wm-alligned_mask.nii.gz
	if [[ ${union_gm_wm} == True ]];then
		3dcalc -a ${alligned_mask} -b ${dir}/${directory}/${subj}_ses-1_${task}_echo-1_part-mag_brain_mask-union.nii.gz  -expr '(a*b)' -prefix ${alligned_mask_union} -overwrite

		3dcalc -a  ${dir}/${directory}/${target}_gm.nii.gz -b ${alligned_mask_union} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_gm_union.nii.gz -overwrite
		3dcalc -a  ${dir}/${directory}/${target}_non-gm.nii.gz -b ${alligned_mask_union} -expr '(a*b)' -prefix ${dir}/${directory}/${target}_non-gm_union.nii.gz -overwrite
	fi
done


if [[ ${residuals} == True ]] && [[ ${union_gm_wm} == True ]];then
	echo "##################################################"
	echo " Using GM in the functional space to mask OC residuals"
	echo "##################################################"
	3dcalc -a ${root}/tedana_${task}_meicaresiduals${volume}/desc-optcom_bold.nii.gz -b ${alligned_mask_union} -expr '(a*b)' -prefix ${dir}/analysis_timeSeries/${target_residuals}_gm_union.nii.gz -overwrite

	3dcalc -a ${root}/tedana_${task}_meicaresiduals${volume}/desc-optcom_bold.nii.gz -b ${non_gm_mask_corrected_union} -expr '(a*b)' -prefix ${dir}/analysis_timeSeries/${target_residuals}_non-gm_union.nii.gz -overwrite

fi

