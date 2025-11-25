#!/bin/bash

module load rocks-freesurfer-5.3.0
module load rocks-afni-latest
module load rocks-fsl-5.0.2.2
module load rocks-mricron-6.2013 

## full/path/to/site/subject_list

subjects_list=subjects_PARK

top_dir=/export/home/ccaballero/ccaballero/data/PARK

align_dir=${top_dir}/review_align_PARK

###############################################################

##--------------------START OF SCRIPT------------------------##

###############################################################

## Get subjects to run

subjects_FOGp_list=subjects_FOGp
subjects_FOGp=$(cat ${subjects_FOGp_list})

subjects_FOGm_list=subjects_FOGm
subjects_FOGm=$(cat ${subjects_FOGm_list})

subjects_CONt_list=subjects_CONt
subjects_CONt=$(cat ${subjects_CONt_list})

#subjects=01C


## move all data to review alignements into the same folder
## SUBJECT LOOP

date

# DELETE DIRECTORY WITH DATASETS
rm -R ${align_dir}
mkdir -p ${align_dir}

cd ${align_dir}

ln -s /opt/afni-2011_12_21_1014/linux_openmp_64/TT_N27+tlrc.BRIK.gz ${align_dir}/TT_N27+tlrc.BRIK.gz
ln -s /opt/afni-2011_12_21_1014/linux_openmp_64/TT_N27+tlrc.HEAD ${align_dir}/TT_N27+tlrc.HEAD

for subject in $subjects_FOGm $subjects_FOGp $subjects_CONt
do

echo "Linking required volumes for ${subject}"
ln -s ${top_dir}/${subject}/proc_RS/T1_brain_unif_al_epi+orig.BRIK.gz ${align_dir}/T1_brain_unif_${subject}_al_epi+orig.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/T1_brain_unif_al_epi+orig.HEAD ${align_dir}/T1_brain_unif_${subject}_al_epi+orig.HEAD

ln -s ${top_dir}/${subject}/proc_RS/mean.pb04.${subject}.tproject+orig.BRIK.gz ${align_dir}/mean.pb04.${subject}.tproject+orig.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/mean.pb04.${subject}.tproject+orig.HEAD ${align_dir}/mean.pb04.${subject}.tproject+orig.HEAD

ln -s ${top_dir}/${subject}/proc_RS/mean.pb04.${subject}.tproject+tlrc.BRIK.gz ${align_dir}/mean.pb04.${subject}.tproject+tlrc.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/mean.pb04.${subject}.tproject+tlrc.HEAD ${align_dir}/mean.pb04.${subject}.tproject+tlrc.HEAD

ln -s ${top_dir}/${subject}/proc_RS/vr_base_min_outlier+orig.BRIK.gz ${align_dir}/${subject}.vr_base_min_outlier+orig.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/vr_base_min_outlier+orig.HEAD ${align_dir}/${subject}.vr_base_min_outlier+orig.HEAD

ln -s ${top_dir}/${subject}/proc_RS/mask_WMe_LowRes+orig.BRIK.gz ${align_dir}/${subject}.mask_WMe_LowRes+orig.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/mask_WMe_LowRes+orig.HEAD ${align_dir}/${subject}.mask_WMe_LowRes+orig.HEAD

ln -s ${top_dir}/${subject}/proc_RS/mask_ventCSF_LowRes+orig.BRIK.gz ${align_dir}/${subject}.mask_ventCSF_LowRes+orig.BRIK.gz
ln -s ${top_dir}/${subject}/proc_RS/mask_ventCSF_LowRes+orig.HEAD ${align_dir}/${subject}.mask_ventCSF_LowRes+orig.HEAD

done

# Open AFNI and avoid the pop-up message for warning of oblique datasets
export AFNI_NO_OBLIQUE_WARNING=YES


cd ${align_dir}

afni -niml -yesplugouts &

# loop across subjects to see the results of

sleep 15


for subject in $subjects_FOGm $subjects_FOGp $subjects_CONt
do

echo "Press ENTER to see ORIG alignements for subject ${subject}:"
read pause

plugout_drive \
-com "SWITCH_UNDERLAY A T1_brain_unif_${subject}_al_epi+orig." \
-com "SWITCH_OVERLAY A mean.pb04.${subject}.tproject+orig" \
-com 'OPEN_WINDOW A.sagittalimage geom=800x600 mont=13x13:1 opacity=4 ifrac=1' \
-com 'OPEN_WINDOW A.axialimage geom=1000x1000 mont=12x12:1 opacity=4 ifrac=1' \
-com 'SET_VIEW A.orig' \
-com 'SET_FUNC_RANGE A.900' \
-com 'SET_PBAR_ALL A.+99 1.0 Spectrum:red_to_blue' \
-com 'SET_PBAR_SIGN A.+' \
-quit

sleep 2

echo "Press ENTER to see WM mask in functional space for subject ${subject}"
read pause


plugout_drive \
-com "SWITCH_UNDERLAY A ${subject}.vr_base_min_outlier+orig." \
-com "SWITCH_OVERLAY A ${subject}.mask_WMe_LowRes+orig." \
-com 'OPEN_WINDOW A.sagittalimage geom=800x600 mont=13x13:1 opacity=4 ifrac=1' \
-com 'OPEN_WINDOW A.axialimage geom=1000x1000 mont=12x12:1 opacity=4 ifrac=1' \
-com 'SET_VIEW A.orig' \
-com 'SET_FUNC_RANGE A.1' \
-com 'SET_PBAR_ALL A.+99 1.0 Spectrum:red_to_blue' \
-com 'SET_PBAR_SIGN A.+' \
-quit

sleep 2

echo "Press ENTER to see ventricle CSF mask in functional space for subject ${subject}"
read pause


plugout_drive \
-com "SWITCH_OVERLAY A ${subject}.mask_ventCSF_LowRes+orig." \
-quit

sleep 2

echo "Press ENTER to see TLRC alignements for subject ${subject}"
read pause


plugout_drive \
-com "SWITCH_UNDERLAY A TT_N27+tlrc" \
-com "SWITCH_OVERLAY A mean.pb04.${subject}.tproject+tlrc" \
-com 'OPEN_WINDOW A.sagittalimage geom=800x600 mont=13x13:1 opacity=4 ifrac=1' \
-com 'OPEN_WINDOW A.axialimage geom=1000x1000 mont=13x13:1 opacity=4 ifrac=1' \
-com 'SET_VIEW A.tlrc' \
-com 'SET_FUNC_RANGE A.900' \
-com 'SET_PBAR_ALL A.+99 1.0 Spectrum:red_to_blue' \
-com 'SET_PBAR_SIGN A.+' \
-quit

sleep 2

echo "Press ENTER to see realignment parameter timeseries and VARS time series for subject ${subject}"
read pause 

1dplot -sepscl -jpeg "motion_${subject}" ${top_dir}/${subject}/proc_RS/motion_demean.1D ${top_dir}/${subject}/proc_RS/motion_${subject}_enorm.1D ${top_dir}/${subject}/proc_RS/motion_${subject}_censor.1D ${top_dir}/${subject}/proc_RS/outcount_${subject}_censor.1D ${top_dir}/${subject}/proc_RS/censor_${subject}_combined_2.1D

n_trs=$(3dinfo -nt mean.pb04.${subject}.tproject+orig)
n_valid_trs=$(1dsum ${top_dir}/${subject}/proc_RS/censor_motion_0.5_${subject}_combined_2.1D) 
echo "The number of valid timepoints is: ${n_valid_trs}"
n_censored_trs=$(echo "${n_trs}-${n_valid_trs}" | bc)
echo "The number of censored timepoints is: ${n_censored_trs}"
# copy all these parameters into a general file
printf "%4s %4s %4s %4s \n" ${subject} ${n_trs} ${n_valid_trs} ${n_censored_trs} >> list_censored_PARK

sleep 2

done

# close AFNI
echo "Hit Enter to Quit AFNI"
read pause
plugout_drive  -com 'QUIT' \
               -quit   

# DELETE DIRECTORY WITH DATASETS
cd ${top_dir}
rm -R ${align_dir}

export AFNI_NO_OBLIQUE_WARNING=NO
