#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.glm.first_level import make_first_level_design_matrix, FirstLevelModel
from nilearn.masking import apply_mask, unmask
from nilearn.image import concat_imgs, index_img, clean_img, get_data
# This script will evaluate the quality of denoising on some task based fmri
# data. The data contains multiple runs, There are several different denoising
# methods - in some case we may want to relate these to several other runs of
# the original (nondenoise data).

# The script should take in
# - a list of files or perhaps a reg exp to match files
# - a list of the timing files associated with the nii data
# - the path the mask we want to use.
# - a flag for the denoising method - method, string
# - a flag to use either 3ddeconvolve or nilearn - use_3dDeconvolve, 1/0
# - a flag for the number of held. out runs - n_hold_out, int
# - a flag for the number of bootstrap iterations - n_iters, int
# - numbner of polynomials to remove - n_poly, int

# NOTE - the DATA should already be scaled so that it is in units of percent signal change.


# Load in the data - these are nii files, so we can load all of them at once and
# they will j ust be stored as references (using nibabel)

source_directory=/
def glm_predict(source_directory,extention,method,TR=2,HRF="spm", n_regressors=6
        )#source directory is a string, extention is file finish pattern,
## source directory mhas to contain .tsv with events time, duration and class
#method is a string


vanilla_files=sorted([os.path.join(root, x) 
    for root,dirs,files in os.walk(source_directory) 
    for x in files if x.endswith(extention)])
vanilla_data_files = []  # List of file paths to the nii data files

denoised_data_files = []  # List of file paths to the denoised nii data files

n_runs = len(vanilla_files)

for i in range(n_runs):
    vanilla_data_files.append(
           apply_mask(vanilla_files[i],
            mask)
           )  # Assuming file_list contains paths to vanilla data files
    #data is already masked
    denoised_data_files.append(
            apply_mask(vanilla_files[i].replace("vanilla", method),
            mask)
            )  # This will load the denoised data that is identical to the og data
    #data is already masked
    # Check that the timing files have n_rows equal to the number of files.
# Loop through the number of permutations (choose n runs to hold out, randomply from whole list, )

n_splits = 100

indices = np.arange(n_runs)
splits = []

n_hold_out = 3
for _ in range(n_splits):
    test_idx = np.random.choice(indices, size=(n_runs//n_hold_out)*2, replace=False)
    train_idx = np.setdiff1d(indices, test_idx)
    splits.append((train_idx, test_idx))

if spc_trans:
    output_vanilla_files=[string.replace("desc","desc-vanillaspc")
            for string in vanilla_files]#creating output nam
    vanilla_data_spc_files=[]
    denoised_data_spc_files=[]
    for i in range(n_runs):
        vanilla_data_spc_files.append(
            (vanilla_data_files[i] - 
            (np.mean(vanilla_data_files[i],axis=0).reshape(-1,1))) /
            (np.mean(vanilla_data_files[i],axis=0).reshape(-1,1))
            )
        denoised_data_spc_files.append(
            (denoised_data_files[i] - 
            (np.mean(denoised_data_files[i],axis=0).reshape(-1,1))) /
            (np.mean(denoised_data_files[i],axis=0).reshape(-1,1))
            )
        if save_spc:
            spc_vanilla = unmask(
            np.transpose(
            vanilla_data_spc_files[i]),
            mask)
            spc_vanilla.to_filename(
            output_vanilla_files[i]
            )
            spc_denoised = unmask(
            np.transpose(
            denoised_data_spc_files[i]),
            mask
            )
            spc_denoised.to_filename(
            output_vanilla_files[i].replace("vanillaspc",method+"spc")
            )
#sometimes I get some voxels with 0, maybe I can add a constant noise
            # Open event tsv
events_timeseries = {}
design_matrix = {}
n_bricks = vanilla_data_files[0].shape[0]
frame_times =np.arange(n_bricks)*TR
events_tsv=sorted([os.path.join(root, x) 
    for root,dirs,files in os.walk(source_directory.split("sub-")[0]) 
    for x in files if x.endswith("events.tsv")])
for i in range(len(events_tsv)):
    tmp = pd.readcsv(
            events_tsv[i],
            sep = '\t',
            header = 0
            )
    events_timeseries[i] = tmp.loc[tmp["trial_type"] != "baseline"]
    design_matrix[i] = make_first_level_design_matrix(frame_times,
            events_timeseries[i],
            drift_model = "polynomial",
            drift_order = "4",
            hrf_model=HRF)
    run = events_tsv[i].split("run-")[-1].split("_")[0]
    design_matrix[i].rename(
        columns={
        "constant":"constant"+run,
        "drift_1":"drift_1"+run,
        "drift_2":"drift_2"+run,
        "drift_3":"drift_3"+run,
        "drift_4":"drift_4"+run
        },inplace=True)
    #as key and the run number as events column name so that then we can acces them and merge them letting columns stay different per run
# splits hold the indicies of the training and testing data for each iteration
# making individual design matrixes

# to create design matrix pd.concat([design_matrix[i],design_matrix[i+n]]) 

#creating confound matrix with polynomials regressor to take out

confound_matrix = make_first_level_design_matrix(
        frame_times,
        drift_model="polynomial",
        drift_order=4,
        hrf_model="spm")

# To avoid redundant polynomial regression, I will do that at the beggining for all
# the loaded volumes, con: more memory, pro: faster 
regressed_confounds_vanilla_files = []


for i in range(n_runs):
#This is messy, but _I am firs detrending the spc timeseries, and then masking it
#This gives me back a matrix time*voxels

    regressed_confounds_vanilla_files.append(
               clean_img(
                   output_vanilla_files[i],
                   confounds=confound_matrix,
                   detrend=False,
                   standardize=False,
                   sample_mask=mask).get_fdata()#detrend polynomials within the brain mask (x*y*z*t)
           )

    regressed_confounds_denoised_files.append(
               clean_img(
                   output_denoised_files[i],
                   confounds=confound_matrix,
                   detrend=False,
                   standardize=False,
                   sample_mask=mask).get_fdata()#detrend polynomials within the brain mask (x*y*z*t)
           )


for test, train in splits:
    # decision point - use nilearn or 3ddeconvolve
    if nilearn:
        
        # taken from https://nilearn.github.io/stable/glm/first_level_model.html
        
        train_vanilla_files = [output_vanilla_files[i] for i in train]

        train_denoised_files = [
                f.replace("vanillaspc",method+"spc") 
                for f in train_vanilla_files
                ]

        test_vanilla_files = [regressed_confounds_vanilla_files[i] for i in test]

        test_denoised_files = [regressed_confounds_denoised_files[i] for i in test]
        
        design_matrices = pd.concat([
            design_matrix[i] for i in train 
            ],ignore_index = True)#Creating one design matrix for all runs

        test_design_matrices = pd.concat([
            design_matrix[i] for i in test
            ],ignore_index = True)# creating design matrix for testing
        
        test_design_matrices = test_design_matrices.values[:,n_regressors]# converting to numpy array for tensor product and taking only the task conditions as regressors
        
        #each onehas individual polynomials, but share the same events
        
        design_matrices = design_matrices.fillna(0)#replace NaN with 0 so that 
        
        # we can run a glm
        
        concatenated_vanilla_images = concat_imgs(train_vanilla_files)
        #nilearn only accepts nilearn objects, so we concatenate in time with the 
        #same order as the design matrices
        
        fmri_train_glm = FirstLevelModel(t_r=TR,
                mask_img=mask,
                standardize=False,
                signal_scaling=False,
                hrf_model=HRF,
                minimize_memory=False)

        fmri_train_glm = fmri_train_glm.fit(concatenated_vanilla_images, design_matrices=design_matrices)#training glm

        regressors = design_matrices.columns.tolist()

        contrast_matrix=np.zeros((n_regressors,len(regressors))

        regressors = regressors[0:n_regressors]# getting only the task regressors
        #which are the first ones

        indices_regressors=list(range(n_regressors))

        for i, idx in enumerate(indices_regressors):
            contrast_matrix[i,idx]=1#contrast matrix to calculate betas

        betas_fmri = get_data(
            fmri_train_glm.compute_contrast(
                contrast_matrix,
                output_type="effect_size"))
             # This is going to return us a pandas dataframe dimentions x,y,z,conditions
        
        print ("Predicting timeseries given betas and test design matrix")

        predicted_timeseries = np.tensordot(
            betas_fmri,
            test_design_matrices,
            axes=([3],[1])
            )#
        test_vanilla_files=np.concatenate(
            test_vanilla_files,
            axis=3
            )

        shape_predicted=predicted_timeseries.shape#saving shape

        shape_test_vanilla=test_vanilla_files.shape#saving shape

        # we will reshape everything to a voxel * Time matrix
        predicted_timeseries = np.reshape(
            predicted_timeseries,
            (np.prod(
                shape_predicted[0:3]),#prod = amount of voxels
                shape_predicted[-1]) # -1 = amount volulmes across time
            )

        test_vanilla_files = np.reshape(
            test_vanilla_files,
            (np.prod(
                shape_test_vanilla[0:3]),#prod = amount of voxels
                shape_test_vanilla[-1]) # -1 = amount volulmes across time
            )
        
        #creating a mask for only voxels of interest




        
        #for regressor in range(n_regressors):
        #    betas_dict[regressors[regressor]] = betas_fmri[regressor,:]#saving each regressor beta in a different key




    # We have the betas - generate the design matrix for the testing data, using the same 3ddecnolve approach

    # Load in the generate design matrix

    # NOTE - You must project out the polynomials from the loaded design matrix (correctly, on a per run basis)
    # before we do the prediction. Why?
    # We have betas that are the results of a model that include polynomials.
    # we are going to have raw data that we have projected the polynomials out of.
    # Therefore, we need to match the design to the data - so we project out the polynomials before we generate the predicted timeseries.
    # I'm 99% sure about this.

    # X: [T, C]      (design matrix, with polys projected out)
    # B: [X, Y, Z, C]  (betas)
    # predicted: [X, Y, Z, T]
    # This creates the prediced timeseries for the entire design matrix
    # as an output that is X, Y, Z, Time
    predicted = np.transpose(np.tensordot(X, B, axes=([1], [3])), (1, 2, 3, 0))

    # actually load in the test data (get_fdata() from nibabel)

    # Project out the polynomials from the data (its like 3dTproject, but you can do it here in code)
    # subselect the polynomial portion of the design matrix
    for run in test:
        # project out polys
        # concatenate the runs together, so that we have a single 4D array

    

   
    # For each voxel, calculate the R2 between the predicted and the actual data

    # Store this in an ndarray, size X, Y, Z, n_iters

# Permutations are done.

# Save out the median of the R2 values. Median is chosen instead of mean, because mean of R2 isn't really valid, I don't think.

# We now have the median R2 for denoised data predicting the held out, non denoised data, for a given number of held out runs.
# If this is too slow, we could mask in 3ddeconvolve or mask/vectorize for nilearn, but have to handle the transforms correctly.
#!/usr/bin/env python3
