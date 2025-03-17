#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=fmriprep_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=09:00:00
#SBATCH --array=0
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/HABLA/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/HABLA/logs/%x_%a_%j.txt


#module load AFNI
module load Apptainer/1.2.4

echo "##################################################"
echo "MEICA JOB ${SLURM_ARRAY_TASK_ID}"
echo "##################################################"


#sub=( "01" "02" "03" "04" "06" "07" "08" "09" "11" "12" "14" "15" "16" "17" "18" "20" "21" )
sub=( "09" )


bids_data=/scratch/mflores/HABLA
out_data=/scratch/mflores/HABLA
bids_derivatives=${out_data}/derivatives/fmriprep_vanilla/
freesurfer_dir=${out_data}/derivatives/fmriprep_freesurfer/sourcedata/freesurfer2/
sub_label=${sub[${SLURM_ARRAY_TASK_ID}]} #Change -N option above when you change sub_label
ses_label=1
#run_label=01
#task_label=rest
sub_outdir=${bids_derivatives}/sub-${sub_label}/
work_dir=/scratch/mflores/temp
#work_dir=${bids_derivatives}/sub-${sub_label}/ses-${ses_label}/work
echo "##################################################"
echo "############# SUBJECT :sub-${sub_label} #############"
echo "##################################################"
mkdir -p $work_dir
mkdir -p $freesurfer_dir
mkdir -p $bids_derivatives
# As a baseline, going to let anything in the work dir stick around because the only reason it would be there is if this script didn't complete normally, in which case, having intermediaries for rerunning might be useful
# rm -r ${work_dir}/fmriprep_23_1_wf/single_subject_${sub_label}_wf/func_preproc_ses_${ses_label}_*_echo_1_wf/

space_avail=`df -h --output=avail /scratch | sed '1d'`
if [ ${space_avail%?} -lt 300 ]
then
    echo "Only $space_avail free space in /scratch, need 300GB to run. Exiting."
    exit
fi


### renaming thermal denoised volumes to follow bids format.

cd ${bids_data}/sub-${sub_label}/ses-${ses_label}/func

echo "cd ${bids_data}/sub-${sub_label}/ses-${ses_label}/func"
for f in *bold.nii.gz; do mv "$f" "$(echo "$f" | sed s/bold/BOLD/)";done

for f in *_vanilla.nii.gz; do mv "$f" "$(echo "$f" | sed s/_vanilla.nii.gz/.nii.gz/)";done


singularity run --cleanenv \
    -B ${bids_data}:/input \
    -B ${bids_derivatives}:/output \
    -B ${work_dir}:/work \
    -B ${freesurfer_dir}:/freesurfer \
    -B /data/mflores/freesurfer/license.txt:/license.txt \
    -B /data/mflores/fmriprep_filters/fmriprep_filter_ses-${ses_label}.json:/fmriprep_filter.json \
    /data/mflores/singularity_images/fmriprep-23.1.3.simg \
    /input /output \
    participant \
    --fs-license-file /license.txt \
    --skip_bids_validation \
    --participant-label $sub_label \
    --bids-filter-file /fmriprep_filter.json \
    --work-dir /work \
    --fs-subjects-dir /freesurfer \
    --output-spaces MNI152NLin6Asym func T1w \
    --nprocs 16 \
    --omp-nthreads 16 \
    --me-output-echos \
    --force-bbr \
    --ignore slicetiming \
    --dummy-scans 10 \
    --fd-spike-threshold 0.3

rm -rf ${work_dir}/fmriprep_23_1_wf/single_subject_${sub_label}_wf/

cd ${bids_data}/sub-${sub_label}/ses-${ses_label}/func


for f in *bold.nii.gz; do mv "$f" "$(echo "$f" | sed s/bold/bold_vanilla/)";done
for f in *BOLD.nii.gz; do mv "$f" "$(echo "$f" | sed s/BOLD/bold/)";done
#@SUMA_Make_Spec_FS -NIFTI -fspath ${freesurfer_dir}sub-${sub_label}/surf/ -sid sub-${sub_label}

