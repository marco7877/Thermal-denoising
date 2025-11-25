#!/bin/bash
######### CVR Images AFNI
# Author:  Cristina Comella
# Version: 1.0
# Date:    13.03.2024
#########

module load afni/latest

#Path/to/subs
echo "path: 
PRJDIR=/bcbl/home/public/CVR/LANGCONN_BIDS"

PRJDIR=/bcbl/home/public/CVR/LANGCONN_BIDS
cvr_mni_folder=${PRJDIR}/images_CVR_MNI



###############################################################

##--------------------START OF SCRIPT------------------------##

###############################################################

#move all data to review alignemets into the same folder 
echo "creating directory cvr mni"
#mkdir -p ${cvr_mni_folder}
#cd ${cvr_mni_folder}


echo "********input_file=${PRJDIR}/langconn_names.txt"
echo "********output_file=${PRJDIR}/output_prueba.txt"


input_file="${PRJDIR}/langconn_names.txt"
output_file="${PRJDIR}/output_prueba.txt"

ls -l ${input_file}

#Read the file line by line

if [[ -L "${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz" ]];then
            ln -s /opt/afni/linux_openmp_64/MNI152_2009_template_SSW.nii.gz ${cvr_mni_folder}/MNI152_2009_template_SSW.nii.gz 
            echo "** ln -s /opt/afni/linux_openmp_64/MNI152_2009_template_SSW.nii.gz ${cvr_mni_folder}/MNI152_2009_template_SSW.nii.gz**"
fi



echo "Reading text file"
while IFS=s';' read -r sub ses class; do
    echo "it enters the loop"

    echo "Sub: $sub" 
    echo "Ses: $ses"
    echo "Classification: $class"

    if [[ "$class" != "good" ]]; then
    echo "***not good the classification"
        continue

    fi
   
   func_preproc_dir=${PRJDIR}/sub-${sub}/ses-${ses}/func_preproc
   echo "Func preproc dir: ${func_preproc_dir}"

   anat_preproc_dir=${PRJDIR}/sub-${sub}/ses-${ses}/anat_preproc
   echo "Anat preproc dir: ${anat_preproc_dir}"

   cvr_results_dir=${func_preproc_dir}/phys2cvr_maps
   echo "CVR dir: ${cvr_results_dir}"
   #echo "$sub $ses $class" >> "$output_file"


    #If link is in the folder eliminate it
    if [[ -L "${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz" ]];then
            rm -f ${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz 
            echo "${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz"
    fi

     
    if [[ -L "${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz" ]];then
            rm -f ${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz
            echo "${cvr_mni_folder}/sub-${sub[i]}_ses-${ses[i]}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz"
    fi


    # if file exist in phys2cvr folder in each subject, copy it in the CVR MNI Folder
   if [[ -e "${cvr_results_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz" ]]; then
        echo "CVR MNI NN is copying into the cvr mni folder" 
        echo "**ln -s ${cvr_results_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz  ${cvr_mni_folder}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz**"
        ln -s ${cvr_results_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz  ${cvr_mni_folder}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz 
   fi
   
    if [[ -e "${anat_preproc_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz" ]]; then
        echo "The Lag MNI is copying into the cvr mni folder"
        echo "**ln -s ${anat_preproc_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz  ${cvr_mni_folder}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz"
        ln -s ${anat_preproc_dir}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz ${cvr_mni_folder}/sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz 
    fi
    


done < "$input_file" 


#2nd Step: Open AFNI and do the automatic part
#Open AFNI and avoid the pop-up for warning of olbique datasets

export AFNI_NO_OBLIQUE_WARNING=YES


cd $cvr_mni_folder


afni -niml -yesplugouts &

sleep 15 

echo "Common parameters for each subject"

read pause


plugout_drive -com "SWITCH_UNDERLAY A MNI152_2009_template_SSW.nii.gz" \
     -com 'OPEN_WINDOW sagittalimage geom=800x600 mont=9x8:2 opacity=7 ifrac=1' \
     -com 'OPEN_WINDOW axialimage geom=1000x1000 mont=9x8:2 opacity=7 ifrac=1' \
     -com 'SET_FUNC_AUTORANGE A.-' \
     -com 'SET_XHAIRS OFF' \
     -com 'SET_PBAR_ALL A.-99 0.30 Spectrum:yellow_to_cyan' \
     -com 'SET_FUNC_RANGE A. 0.30' \
     -quit

input_file="${PRJDIR}/langconn_names.txt"
ls -l ${input_file}
echo "Reading text file"

echo "Reading text file"
while IFS=s';' read -r sub ses class; do
    echo "it enters the CVR loop"

    echo "Sub: $sub" 
    echo "Ses: $ses"
    echo "Classification: $class"

    if [[ "$class" != "good" ]]; then
    echo "***not good the classification"
        continue

    fi
     
     # Espera unos segundos para asegurarte de que las ventanas se han abierto completamente

     plugout_drive -com "SWITCH_OVERLAY A sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_MNI_NN.nii.gz" \
               -com "SAVE_JPEG axialimage sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_NN_axial.jpg" \
               -com "SAVE_JPEG sagittalimage sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_NN_sagittal.jpg "\
               -quit

    sleep 8

done < "$input_file" 


plugout_drive -com 'SWITCH_UNDERLAY A MNI152_2009_template_SSW.nii.gz' \
     -com 'OPEN_WINDOW sagittalimage geom=800x600 mont=9x8:2 opacity=7 ifrac=1' \
     -com 'OPEN_WINDOW axialimage geom=1000x1000 mont=9x8:2 opacity=7 ifrac=1' \
     -com 'SET_FUNC_AUTORANGE A.-' \
     -com 'SET_XHAIRS OFF' \
     -com 'SET_PBAR_ALL A.-99 5.0 Viridis'  \
     -com 'SET_FUNC_RANGE A. 5' \
     -quit



echo "Reading text file"
while IFS=s';' read -r sub ses class; do
    echo "it enters the Lag loop"

    echo "Sub: $sub" 
    echo "Ses: $ses"
    echo "Classification: $class"

    if [[ "$class" != "good" ]]; then
    echo "***not good the classification"
        continue

    fi



     plugout_drive -com "SWITCH_OVERLAY A sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_MNI_NN.nii.gz" \
               -com 'SET_PBAR_ALL A.-99 5.0 Viridis'  \
               -com "SAVE_JPEG axialimage sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_NN_axial.jpg" \
               -com "SAVE_JPEG sagittalimage sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_NN_sagittal.jpg" \
               -quit

    sleep 8

done < "$input_file" 


plugout_drive -quit -com 'QUIT' \

##Step 3: move images

while IFS=s';' read -r sub ses class; do
    echo "Moving all CVR files and Lag files in eahc folder"

    echo "Sub: $sub" 
    echo "Ses: $ses"
    echo "Classification: $class"

    if [[ "$class" != "good" ]]; then
    echo "***not good the classification"
        continue

    fi
    
    mv sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_NN_axial.jpg ${cvr_mni_folder}/CVR_axial2 
    mv sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_cvr_NN_sagittal.jpg ${cvr_mni_folder}/CVR_sagittal2 
    mv sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_NN_axial.jpg ${cvr_mni_folder}/LAG_axial2 
    mv sub-${sub}_ses-${ses}_task-BH_run-1_optcom_bold_sm_lag_NN_sagittal.jpg ${cvr_mni_folder}/LAG_sagittal2
    

done < "$input_file" 