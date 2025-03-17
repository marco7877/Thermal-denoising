#!/bin/bash -l
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=GLM_vanilla_habla_OC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=30G
#SBATCH --array=0-26
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mflores@bcbl.eu
#SBATCH --error=/scratch/mflores/HABLA/logs/%x_%a_%j.txt
#SBATCH --output=/scratch/mflores/HABLA/logs/%x_%a_%j.txt

module load AFNI/23.0.03-foss-2022a

${SLURM_ARRAY_TASK_ID}
file_condition=/scratch/mflores/HABLA/sequences.csv
counter=0
while read line;do
	LINE=($line)
	sub=sub-${LINE[0]}
	task=task-${LINE[1]}
	seqs=( ${LINE[2]} ${LINE[3]} )
	((counter++))
	if [[ $counter -eq ${SLURM_ARRAY_TASK_ID} ]];then
		echo "Reached permutatuion target value"
		break 1
	fi
done < $file_condition

echo"########################################"
echo" subject $sub"
echo" task $task"
echo" sequences $seqs "
echo"########################################"
#sub=sub-12
#task=task-SENT
runs=( "1" "2" )
#seqs=( "1" "2" )
tree=meica

dir="/scratch/mflores/HABLA/derivatives/fmriprep_vanilla"
out_dir="/scratch/mflores/HABLA/derivatives/fmriprep_vanilla"
MultiEcho=True
denoised=False
doMot=True
#MultiEcho=${1:-True}
#doMot=${2:-True}
#START
fwhm=True
suffix=${task}
mot_threshold=.3
echo "${runs[@]}"
if [[ "${MultiEcho}" == True && "${denoised}" == True ]];then
	for run in "${runs[@]}";do
		MEOC=tedana_${task}_run-${run}_${tree}
		echo "******************************"
		echo "WARNING: using brain mask from run 1!!"
		an_dir="${dir}/${sub}/ses-1/func"
		mask=${an_dir}/${sub}_ses-1_${task}_run-1_part-mag_space-MNI152NLin6Asym_desc-brain_mask.nii.gz
		epi_input=${an_dir}/${MEOC}/desc-denoised_bold
		#epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold
		if [[ ! -e "${epi_input}.nii.gz" ]];then
			echo "******************************"
			echo "WARNING: EPI FILE DOES NOT EXIST!"
			echo "${epi_input}"
			echo "******************************"
			exit 1
		fi
		target=${an_dir}/${subj}_glm_${task}_OC
		mkdir -p ${target}
	done
fi

if [[ "${MultiEcho}" == True && "${denoised}" == False ]];then
	for run in "${runs[@]}";do
		echo "******************************"
		echo "WARNING: using brain mask from run 1!!"
		an_dir="${dir}/${sub}/ses-1/func"
		out_an_dir="${out_dir}/${sub}/ses-1/func"
		mask=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-brain_mask.nii.gz
		epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-preproc_bold
		if [[ ! -e "${epi_input}.nii.gz" ]];then
			echo "******************************"
			echo "WARNING: EPI FILE DOES NOT EXIST!"
			echo "${epi_input}"
			echo "******************************"
			exit 1
		fi
		target=${out_an_dir}/${sub}_glm_${task}_OC_${mot_threshold}
		mkdir -p ${target}
	done
fi


if [[ "${MultiEcho}" == False ]]; then
	echo "Analizing Single-echo data:"
	for run in "${runs[@]}";do
		an_dir="${dir}/${sub}/ses-1/func"
		out_an_dir="${out_dir}/${sub}/ses-1/func"
		mask=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-brain_mask.nii.gz
		epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_echo-2_part-mag_desc-preproc_bold
		if [[ ! -e "${epi_input}.nii.gz" ]];then
			echo "******************************"
			echo "WARNING: EPI FILE DOES NOT EXIST!"
			echo "${epi_input}"
			echo "******************************"
			exit 1
		fi
	done
	target=${out_an_dir}/${sub}_glm_${task}_SE_${mot_threshold}
	mkdir -p ${target}
fi

if [[ ! -e ${mask} ]];then
	echo "******************************"
	echo "WARNING: MASK FILE DOES NOT EXIST!"
	echo "${mask}"
	echo "******************************"
	exit 1
fi


glm_input1D=${out_an_dir}/${sub}_ses-1_${task}_runs_part-mag_desc-confounds_timeseries
if [[ "${doMot}" == True ]];then
	echo "Considering motion parameters"
	if [[ ! -e ${glm_input1D}.1D ]];then
		for run in "${runs[@]}";do
			motFile=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-confounds_timeseries.tsv
			File1D=${out_an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-confounds_timeseries
			#############################
			#TODO: Write a function to do this and extendit to create a 1.D for n runs
			#############################
			headers_to_find=("rot_z" "rot_x" "rot_y" "trans_z" "trans_x" "trans_y" "rot_z_derivative1" "rot_x_derivative1" "rot_y_derivative1" "trans_z_derivative1" "trans_x_derivative1" "trans_y_derivative1")

			echo "-----------------------------------"
			echo "Extracting roll, pitch, yaw, dS, dL, dP from fmriprep :${motFile}"
			echo "where the corresponding column names are: ${headers_to_find}"
			echo "-----------------------------------"

			# Read the first line of the file (header line)
			header_line=$(head -n 1 "$motFile") 
			# Split the header line into an array using tab as delimiter
			IFS=$'\t' read -r -a header_array <<< "$header_line"
			indexes=()
			# Loop through headers_to_find and get their indexes
			for header in "${headers_to_find[@]}";do
				for i in "${!header_array[@]}";do
					if [[ "${header_array[$i]}" == "$header" ]];then
						echo "Header: $header is at index: $i"
						indexes+=($i)
					fi
				done
			done

			echo "${indexes[@]}" 
			awk -v idxs="${indexes[*]}" '
			BEGIN { split(idxs, idx_arr, " "); }
			NR > 1 {
			for(i in idx_arr) {
				printf "%s ", $(idx_arr[i]+1);
			}
			print "";
		}' FS='\t' "$motFile" | sed 's/n\/a/0/g' > "${File1D}.1D"
		cat ${File1D}.1D >> ${glm_input1D}.1D

	done 
fi
fi
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
if [[ ! -e ${glm_input1D}.1D ]];then
	echo "******************************"
	echo "WARNING: MOTION CONFOUNDS FILE DOES NOT EXIST!"
	echo "${glm_input1D}"
	echo "******************************"
	exit 1

fi
if [[ "${doMot}" == True ]];then
	if [[ ! -e ${glm_input1D}_Demean.1D ]];then
		echo "Demean 1D file"
		1d_tool.py -infile ${glm_input1D}.1D -demean -write ${glm_input1D}_Demean.1D
	fi
fi

if [[ "${mot_threshold}" != 0 ]];then
	echo "******************************"
	echo "WARNING: Creating motion based censoring"
	echo "******************************"
	der_mot="${target}/${sub}_ses-1_${task}_runs_part-mag_desc-derivatives"
	echo "${der_mot}"
	cut -d " " -f7-12 ${glm_input1D}.1D > ${der_mot}_timeseries.1D
	echo "L1 censoring frames based on motion parameters (derivatives)"
	echo "Derivatives stored as: ${der_mot}_timeseries.1D " 
	echo "Volumes to be censored: >= ${mot_threshold} mm"
	# TODO renombrar el archiv ofinal y agregar el nivel de censor. modularizar tambiééen 
	# para poder tener estadisticas de 0.3 y 0.5 

	1d_tool.py -infile ${der_mot}_timeseries.1D -set_nruns 2 \
		-censor_prev_TR -censor_next_TR \
		-censor_motion ${mot_threshold} ${der_mot}\
		-show_censor_count \
		-overwrite
		for seq in "${seqs[@]}";do
			cat /scratch/mflores/HABLA/1D_HABLA/${seq}_${task}_multivariate.1D >> /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_onsets_${mot_threshold}.1D
		done

fi

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------


if [[ "${MultiEcho}" == True && "${denoised}" == True ]];then

	for run in "${runs[@]}";do
		MEOC=tedana_${task}_run-${run}_${tree}
		echo "******************************"
		echo "WARNING: moving functional space to T1, to MNI"
		ref_MNI=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_boldref.nii.gz
		epi_input=${an_dir}/${MEOC}/desc-denoised_bold.nii.gz
		denoised_mni=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_denoised_bold.nii.gz
		EPI_2_T1=${an_dir}/${sub}_ses-1_${task}_run-${run}_from-T1w_to-scanner_mode-image_xfm.txt
		T1_2_MNI=${dir}/${sub}/ses-1/anat/${sub}_ses-1_acq-uniclean_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
		echo "Running ANTS"

		echo "antsApplyTransforms -i ${epi_input} -o ${denoised_mni} -r ${ref_MNI} -t ${EPI_2_T1} -t ${T1_2_MNI} -n BSpline[5] -e 3 -v 1"

		antsApplyTransforms -i ${epi_input} -o ${denoised_mni} -r ${ref_MNI} -t ${EPI_2_T1} -t ${T1_2_MNI} -n BSpline[5] -e 3 -v 1

		echo "------------------------------ Changing header"
		echo " Template Space: ORIG -> MNI"
		echo " Dataset Type: Echo PLanar -> Anat Bucket"
		3drefit -space MNI -'abuc' ${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_denoised_bold.nii.gz
	done
fi


#-----------------------------------------------
#TODO: not binary fut gausian diameter

if [[ "${fwhm}" == True && "${denoised}" == True ]];then
	for run in "${runs[@]}";do
		epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_denoised_bold
		if [[ ! -e ${epi_input}_fwhm.nii.gz ]];then
			echo "Smoothing EPI file"
			3dBlurInMask -input ${epi_input}.nii.gz -FWHM 5 -mask ${mask} -prefix ${epi_input}_fwhm.nii.gz -overwrite
		fi
	done
fi
if [[ "${fwhm}" == True && "${denoised}" == False ]];then
	for run in "${runs[@]}";do
		epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-preproc_bold
		if [[ ! -e ${epi_input}_fwhm.nii.gz ]];then
			echo "Smoothing EPI file"
			3dBlurInMask -input ${epi_input}.nii.gz -FWHM 5 -mask ${mask} -prefix ${epi_input}_fwhm.nii.gz -overwrite
		fi
	done
fi


echo "-----------------------------------"
echo "Demeaning"
glm_input=()
if [[ "${denoised}" == True ]];then
for run in "${runs[@]}";do
	epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_space-MNI152NLin6Asym_denoised_bold_fwhm
	3dTstat -mean -prefix ${epi_input}_mean.nii.gz ${epi_input}.nii.gz -overwrite
	echo "Calculating signal percent change"
	3dcalc -a ${epi_input}.nii.gz\
		-b ${epi_input}_mean.nii.gz \
		-c ${mask} \
		-expr 'c*min(200, a/b*100)*step(a)*step(b)' \
		-prefix ${epi_input}_spc.nii.gz \
		-overwrite 
	glm_input=${glm_input}" "${epi_input}_spc.nii.gz
done
else
for run in "${runs[@]}";do
	epi_input=${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_desc-preproc_bold_fwhm
	3dTstat -mean -prefix ${epi_input}_mean.nii.gz ${epi_input}.nii.gz -overwrite
	echo "Calculating signal percent change"
	3dcalc -a ${epi_input}.nii.gz\
		-b ${epi_input}_mean.nii.gz \
		-c ${mask} \
		-expr 'c*min(200, a/b*100)*step(a)*step(b)' \
		-prefix ${epi_input}_spc.nii.gz \
		-overwrite 
	glm_input=${glm_input}" "${epi_input}_spc.nii.gz
done
fi



#-----------------------------------------------------------------------------------------
# Prepare task onset times 
for seq in "${seqs[@]}";do
	cat /scratch/mflores/HABLA/1D_HABLA/${seq}_${task}_control.1D >> /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_control.1D
	cat /scratch/mflores/HABLA/1D_HABLA/${seq}_${task}_fricative.1D >> /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_fricative.1D
	cat /scratch/mflores/HABLA/1D_HABLA/${seq}_${task}_plosive.1D >> /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_plosive.1D
	cat /scratch/mflores/HABLA/1D_HABLA/${seq}_${task}_nasal.1D >> /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_nasal.1D
done


#-----------------------------------------------------------------------------------------
# Prepare 3dDeconvolve file

AFNI_GLTSYM_PRINT=YES

echo -e "\033[0;32m++ STEP (1) Compute statistical maps of activation per event type \033[0m"
echo -e "\033[0;32m++ ===============================================================\033[0m"
run3dDeconvolve="3dDeconvolve -overwrite -input${glm_input} -mask ${mask}"
if [[ ${mot_threshold} != 0 ]];then
	run3dDeconvolve="${run3dDeconvolve}  -censor ${der_mot}_censor.1D"
run3dDeconvolve="${run3dDeconvolve} -polort A"
fi

# We have 4 levels: control, fricative, plosive, nasal
#
#We have 3 conditions: sylable, word, sentence
#
#We are taking always 

#declare -A onsetNames
#onsetNames=([task-SYLAB]=0 [task-WORD]=1 [task-SENT]=2)
run3dDeconvolve="${run3dDeconvolve} -num_stimts 4"

run3dDeconvolve="${run3dDeconvolve} -stim_label 1 CONTR -stim_times 1 /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_control.1D 'SPMG1'"
run3dDeconvolve="${run3dDeconvolve} -stim_label 2 FRIC -stim_times 2 /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_fricative.1D 'SPMG1'"
# Using default
run3dDeconvolve="${run3dDeconvolve} -stim_label 3 PLOS -stim_times 3 /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_plosive.1D 'SPMG1'"
run3dDeconvolve="${run3dDeconvolve} -stim_label 4 NAS -stim_times 4 /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_nasal.1D 'SPMG1'"
run3dDeconvolve="${run3dDeconvolve} -ortvec ${glm_input1D}_Demean.1D motdemean"
run3dDeconvolve="${run3dDeconvolve} -num_glt 10"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC -CONTR' -glt_label 1 FRIC_vs_CONTR"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +PLOS -CONTR' -glt_label 2 PLOS_vs_CONTR"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +NAS -CONTR' -glt_label 3 NAS_vs_CONTR"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC -PLOS' -glt_label 4 FRIC_vs_PLOS"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC -NAS' -glt_label 5 FRIC_vs_NAS"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +PLOS -NAS' -glt_label 6 PLOS_vs_NAS"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC +PLOS +NAS' -glt_label 7 ALL_T"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC \ +PLOS \ +NAS' -glt_label 8 ALL_F"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC +PLOS +NAS - 3*CONTR' -glt_label 9 ALL_vs_CONTR_T"
run3dDeconvolve="${run3dDeconvolve} -gltsym 'SYM: +FRIC -CONTR \ +PLOS -CONTR \ +NAS -CONTR' -glt_label 10 ALL_vs_CONTR_F"


# run3dDeconvolve="${run3dDeconvolve} -ortvec ${prepTaskDir}/motDeriv.1D motderiv"   

run3dDeconvolve="${run3dDeconvolve} -tout -fout"
run3dDeconvolve="${run3dDeconvolve} -bucket ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.OLS.nii.gz"
run3dDeconvolve="${run3dDeconvolve} -x1D ${target}/${sub}_ses-1_${task}_run-${run}_part-mag__bold.spm.xmat.1D"
run3dDeconvolve="${run3dDeconvolve} -xjpeg ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.xmat.jpg"
run3dDeconvolve="${run3dDeconvolve} -fitts ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.fitts.nii.gz"
run3dDeconvolve="${run3dDeconvolve} -errts ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.errts.nii.gz"


#Run 3dDeconvolve
echo "${run3dDeconvolve}"
eval ${run3dDeconvolve}

# Run 3dREMLfit with the same model to account for serial correlations
# ==================================================================
echo -e "\033[0;32m++ STEP (2) Run 3dREMLfit with same model \033[0m"
echo -e "\033[0;32m++ =======================================\033[0m"
3dREMLfit -overwrite -matrix ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.xmat.1D \
	-mask "${mask}"         \
	-input "${glm_input}"   \
	-tout -fout -rout -verb   \
	-Rbuck  ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.REML.nii.gz \
	-Rfitts  ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.Rfitts.nii.gz \
	-Rerrts  ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.Rerrts.nii.gz


echo "*******************************************************"
echo "*******************************************************"
echo "** GLM analysis of ${sub} ${task} COMPLETE!  **"
echo "*******************************************************"
echo "*******************************************************"

3drefit -space MNI -'abuc' ${target}/${sub}_ses-1_${task}_run-${run}_part-mag__bold.spm.REML.nii.gz
3drefit -space MNI -'abuc' ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.Rfitts.nii.gz
3drefit -space MNI -'abuc' ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.Rerrts.nii.gz
3drefit -space MNI -'abuc' ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_bold.spm.OLS.nii.gz

echo "making symbolic link for reference bold image"
echo "making symbolic link for reference bold image"

for run in "${runs[@]}";do
	cp -s ${an_dir}/${sub}_ses-1_${task}_run-${run}_part-mag_boldref.nii.gz ${target}/${sub}_ses-1_${task}_run-${run}_part-mag_boldref.nii.gz
	done
echo "----------------------------------------"
echo "Removing temporal files:"
echo "/scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_control.1D"
echo "/scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_fricative.1D"
echo "/scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_plosive.1D"
echo "/scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_nasal.1D"

rm /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_control.1D
rm /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_fricative.1D
rm /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_plosive.1D
rm /scratch/mflores/HABLA/1D_HABLA/${sub}_${task}_nasal.1D



