addpath("/home/mflores/mppca_denoise/");
sub="SUBJECT";
disp(append("Analizing subject: ",sub));
file=append("/scratch/mflores/Resting_State/",sub,"/ses-1/func_preproc_mppca/",sub,"_ses-1_TASK_echo-");
part="_part-mag_bold_mppca_dsd";
ext="HALF.nii.gz";
mag1=append(file,"1",part,ext);
mag2=append(file,"2",part,ext);
mag3=append(file,"3",part,ext);
mag4=append(file,"4",part,ext);
part_out="_part-mag_bold_mppca";
% loading nifti files 
disp(append("reading NIFTI: ",mag1));
nii1=abs(single(niftiread(mag1)));
disp(append("reading NIFTI: ",mag2));
nii2=abs(single(niftiread(mag2)));
disp(append("reading NIFTI: ",mag3));
nii3=abs(single(niftiread(mag3)));
disp(append("reading NIFTI: ",mag4));
nii4=abs(single(niftiread(mag4)));
% loading nifti metadata
nii1_info=niftiinfo(mag1);
nii2_info=niftiinfo(mag2);
nii3_info=niftiinfo(mag3);
nii4_info=niftiinfo(mag4);
% denoising data 
[nii1_denoised,Sigma1]=MPdenoising(nii1,[],[4 4 4],["full"])
[nii2_denoised,Sigma2]=MPdenoising(nii2,[],[4 4 4],["full"])
[nii3_denoised,Sigma3]=MPdenoising(nii3,[],[4 4 4],["full"])
[nii4_denoised,Sigma4]=MPdenoising(nii4,[],[4 4 4],["full"])
% saving denoised files
disp(append("Saving denoised matrix as NIFTI: ",file,"1",part_out,ext));
niftiwrite(nii1_denoised,append(file,"1",part_out,ext),nii1_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"2",part_out,ext));
niftiwrite(nii2_denoised,append(file,"2",part_out,ext),nii2_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"3",part_out,ext));
niftiwrite(nii3_denoised,append(file,"3",part_out,ext),nii3_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"4",part_out,ext));
niftiwrite(nii4_denoised,append(file,"4",part_out,ext),nii4_info);
