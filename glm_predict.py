#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import os
from nibabel import Nifti1Image
from nilearn.image import load_img, resample_to_img
from nilearn.masking import apply_mask, unmask
from nilearn.plotting import plot_epi, plot_stat_map, show
from scipy.stats import pearsonr

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





 sorted([os.path.join(root, x) 
    for root,dirs,files in os.walk(bids_dir) 
    for x in files if x.endswith("echo-1_part-mag_sbref.nii.gz")])

def glm_predict(source_directory,extention,method
        )#source directory is a string, extention is file finish pattern,
#method is a string
vanilla_files=sorted([os.path.join(root, x) 
    for root,dirs,files in os.walk(source_directory) 
    for x in files if x.endswith(extention)])
vanilla_data_files = []  # List of file paths to the nii data files
denoised_data_files = []  # List of file paths to the denoised nii data files
for file in vanilla_files:
    vanilla_data_files.append(
        nib.load(file)
    )  # Assuming file_list contains paths to vanilla data files
    denoised_data_files.append(
        nib.load(file.replace("vanilla", method))
    )  # This will load the denoised data that is identical to the og data

# Check that the timing files have n_rows equal to the number of files.


# Loop through the number of permutations (choose n runs to hold out, randomply from whole list, )
import numpy as np

n_runs = len(vanilla_data_files)
n_splits = 100

indices = np.arange(n_runs)
splits = []

for _ in range(n_iters):
    test_idx = np.random.choice(indices, size=n_hold_out, replace=False)
    train_idx = np.setdiff1d(indices, test_idx)
    splits.append((train_idx, test_idx))

# splits hold the indicies of the training and testing data for each iteration
for train, test in splits:
    # decision point - use nilearn or 3ddeconvolve
    if use_3dDeconvolve:
        # Use 3ddeconvolve
        # the model is simple - just a design matrix with conditions and (per run) polynomials.
        # This will be a system call to 3ddeconvolve, so we need to set up the command line arguments
        # and then call it using subprocess or os.system. use n_poly to set the number of polynomials to remove.
        os.subprocess()

        # load in the betas
        nib.load().get_fdata()  # This will load the betas from the 3ddeconvolve output

        # make sure you know what the shape is, becaeuse we need that to make the predicted timeseries.

        # Done, final betas here represent the (percent dignal change) calulated from how every many runs remained after the held out runs.

    else:
        # Use nilearn (we can also use direct linear algebra if we want?)
        # There will still be a 3ddeconovle run to set up the design matrix

        os.subprocess()

        # Load in the designa matrix.
        # it is special text, but should be reasonble to fgit.

        # fit the glm
        # taken from https://nilearn.github.io/stable/glm/first_level_model.html
        # I have no idea how to do a glm in nilearn, but someone you know does it seems.
        from nilearn.glm.first_level import FirstLevelModel

        fmri_glm = FirstLevelModel()
        fmri_glm = fmri_glm.fit(subject_data, design_matrices=design_matrices)

        # get the betas on the testing data
        # its somewhere in fmri_glm.

        # Done, final betas here represent the (percent dignal change) calulated from how every many runs remained after the held out runs.

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
