#!/bin/bash
#$ -m be
#$ -m be
#$ -M m.flores@bcbl.eu
#$ -S /bin/bash
module load afni/stable

method=$1
sub=$2
tsk=$3
analysis_windows=analysis_03-24_
repo="/bcbl/home/public/MarcoMotion/scripts/HABLA_SPiN/"
target="/bcbl/home/public/MarcoMotion/Habla_restingState/${analysis_windows}/"
preproc="func_preproc_${method}/"


mkdir -p ${target}


origin="/bcbl/home/public/MarcoMotion/Habla_restingState/sub-00${sub}/ses-1/func_preproc_vanilla/"
root="/bcbl/home/public/MarcoMotion/Habla_restingState/sub-00${sub}/ses-1/${preproc}"

echo analizing subject sub-00${sub}
echo analizing task ${tsk}

mask=${origin}sub-00${sub}_ses-1_${tsk}_echo-1_part-mag_brain_mask.nii.gz
ln_mask=${target}sub-00${sub}_ses-1_${tsk}_echo-1_part-mag_brain_mask.nii.gz

ln -s ${mask} ${ln_mask}

polort=5 # Legendre polynomials for detrending
#nhbd_size=5 # size of neighbourhood in mm for computing smoothness with 3dLocalstat


echo "##############################"
echo "tSNR calculation"
echo "##############################"

echo "Compute voxelwise mean"
echo "3dTstat -mean -prefix  ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mean.nii.gz \
		${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mcf_al.nii.gz -overwrite"

3dTstat -mean -prefix  ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mean.nii.gz \
		${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mcf_al.nii.gz -overwrite
echo "Detrend dataset"
echo "3dTproject -polort ${polort} -prefix ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_dt.nii.gz \
		${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mcf_al.nii.gz -overwrite"

3dTproject -polort ${polort} -prefix ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_dt.nii.gz -input  ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mcf_al.nii.gz -overwrite
echo "Put back the mean"

3dcalc -a ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_dt.nii.gz \
		-b ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mean.nii.gz \
		-expr 'a+b' -prefix ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_dt.nii.gz -overwrite

echo "Compute voxelwise standard deviation"
3dTstat -stdevNOD -prefix  ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_std.nii.gz \
		${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_dt.nii.gz -overwrite

echo "Compute tSNR map"
3dcalc -a ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_mean.nii.gz \
		-b ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_std.nii.gz \
		-m ${mask} \
		-expr 'm*(a/b)' \
		-prefix ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_tsnr.nii.gz \
		-overwrite

echo "##############################"
echo "tSNR symbolic links"
echo "##############################"
ln -s ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_tsnr.nii.gz ${target}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${method}_tsnr.nii.gz

#Estimating spatial smoothness per voxel area = to a Sphere with 2x voxel size
#command 3dLocalStat -nbhd 'SPHERE(voxelsize)' -stat FWHM -mask ${mask} -prefix -overwrite

echo "#############################################################################"
echo " Computting global smoothness with 3dFWHMx on raw data and after realignment "
echo "#############################################################################"

for ext_data in ${method} ${method}_mcf_al
do
3dFWHMx  -input ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}.nii.gz \
		-mask ${mask} -out ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}.fwhm -overwrite

done

echo "###########################################################"
echo " Computting local smoothness maps with 3dLocalstat FWHMbar "
echo "###########################################################"

for ext_data in ${method} ${method}_mcf_al
do

nhbd_size=$( 3dinfo -adi ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}.nii.gz )

	3dLocalstat -nbhd "RHDD(5)" -stat FWHMbar -mask ${mask} \
		-prefix  ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhm.nii.gz \
		${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}.nii.gz -overwrite

		ln -s ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhm.nii.gz ${target}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhm.nii.gz

		3dTstat -median -prefix ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhmstats.nii.gz \
			${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhm.nii.gz -overwrite

			ln -s ${root}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhmstats.nii.gz ${target}sub-00${sub}_ses-1_${tsk}_echo-2_part-mag_bold_${ext_data}_fwhmstats.nii.gz


done

