#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=ReconAll2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=11:00:00
#SBATCH --mem=80G
#SBATCH --array=1-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/Resting_State/logs/%x_%a_%j.txt

module load AFNI
module load FreeSurfer

echo "##################################################"
echo "Reconall FRESURFER JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"
list_subj=( "002" "003" )
counter=0
for n_sub in "${list_subj[@]}" 
do
			((counter++))
			if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]]; then
				echo "Reached expected permutation value: ${SLURM_ARRAY_TASK_ID}"
				break 2

			fi
done



echo "##################################################"
echo "############# SUBJECT :${n_sub}      #############"
echo "##################################################"


PRJDIR=/scratch/mflores/Resting_State

subj=sub-${n_sub}

# run Freesurfer recon-all to separate white matter from gray one. For this we are using both the T1 and the T2
# TODO study Freesurfer
recon-all -i ${PRJDIR}/${subj}/ses-1/anat/${subj}_ses-1_acq-uniclean_T1w.nii.gz -T2 ${PRJDIR}/${subj}/ses-1/anat/${subj}_ses-1_T2w.nii.gz -3T -all -s ${subj} -sd ${PRJDIR}/freesurfer
# Import Freesurfer results into SUMA-land 
# TODO study SUMA
@SUMA_Make_Spec_FS -fs_setup -NIFTI -sid ${subj} -fspath ${PRJDIR}/freesurfer/${subj}

#Computing GRAY MATTER and SUBCORTICAL masks from SUMA + FS

labeled_T1=${PRJDIR}/freesurfer/${subj}/SUMA/aparc.a2009s+aseg.nii.gz

gm_mask=${subj}_ses-1_gm-mask_T1w.nii.gz
slf_mask=${subj}_ses-1_slf-mask_T1w.nii.gz
slb_mask=${subj}_ses-1_slb-mask_T1w.nii.gz
sl_mask=${subj}_ses-1_sl-mask_T1w.nii.gz
srf_mask=${subj}_ses-1_srf-mask_T1w.nii.gz
srb_mask=${subj}_ses-1_srb-mask_T1w.nii.gz
sr_mask=${subj}_ses-1_sr-mask_T1w.nii.gz
sc_mask=${subj}_ses-1_sc-mask_T1w.nii.gz

# gm
3dcalc -a lh.ribbon.nii.gz -b rh.ribbon.nii.gz -expr '(a+b)' -prefix ${gm_mask} -overwrite

# subcortical left hemisphere mask

3dcalc -a $labeled_T1 -expr 'within(a,9,13)' -prefix ${slf_mask} -overwrite

3dcalc -a $labeled_T1 -expr 'within(a,17,19)' -prefix ${slb_mask} -overwrite 

3dcalc -a ${slf_mask} -b ${slb_mask} -expr '(a+b)' -prefix ${slf_mask} -overwrite

rm ${slf_mask} ${slb_mask}

# subcortical right hemisphere mask

3dcalc -a $labeled_T1 -expr 'within(a,48,51)' -prefix ${srf_mask} -overwrite

3dcalc -a $labeled_T1 -expr 'within(a,51,55)' -prefix ${srb_mask} -overwrite 

3dcalc -a ${srf_mask} -b ${srb_mask} -expr '(a+b)' -prefix ${sr_mask} -overwrite

rm ${srf_mask} ${srb_mask}

# subcortical mask

3dcalc -a ${sr_mask} -b ${sl_mask} -expr '(a+b)' -prefix ${sc_mask} -overwrite

# target mask

3dcalc -a ${gm_mask} -b ${sc_mask} -expr '(a+b)' -prefix ${subj}_ses-1_analysis-mask_T1w.nii.gz



