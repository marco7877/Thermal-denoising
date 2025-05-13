addpath("/home/mflores/Tensor-MP-PCA/");
sub="SUBJECT";
disp(append("Analizing subject: ",sub));

file=append("/scratch/mflores/Resting_State/",sub,"/ses-1/func_preproc_tmmpca/",sub,"_ses-1_TASK_echo-");
part="_part-mag_bold_tmmpca_dsd";
ext="HALF.nii.gz";
mag1=append(file,"1",part,ext);
mag2=append(file,"2",part,ext);
mag3=append(file,"3",part,ext);
mag4=append(file,"4",part,ext);

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
%creating tensor
disp(append("Creating tensor from magnitude parts"));
tensor_nii=cat(5,nii1,nii2,nii3,nii4);
%running tensor MP-PCA
disp(append("Denoising magnitude valumes from subject: ",sub));
[denoised,Sigma2,P,SNR_gain]=denoise_recursive_tensor(tensor_nii, [4 4], indices={1:3 4 5});
% deconcatenating nifti files
disp("Denoising finishes");
nii1_denoised=denoised(:,:,:,:,1);
nii2_denoised=denoised(:,:,:,:,2);
nii3_denoised=denoised(:,:,:,:,3);
nii4_denoised=denoised(:,:,:,:,4);
%saving denoised parts
%file=strrep(file,"_dsd","");
part_out="_part-mag_bold_tmmpca";
disp(append("Saving denoised matrix as NIFTI: ",file,"1",part_out,ext));
niftiwrite(nii1_denoised,append(file,"1",part_out,ext),nii1_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"2",part_out,ext));
niftiwrite(nii2_denoised,append(file,"2",part_out,ext),nii2_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"3",part_out,ext));
niftiwrite(nii3_denoised,append(file,"3",part_out,ext),nii3_info);
disp(append("Saving denoised matrix as NIFTI: ",file,"4",part_out,ext));
niftiwrite(nii4_denoised,append(file,"4",part_out,ext),nii4_info);
