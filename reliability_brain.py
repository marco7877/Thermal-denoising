#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""
import numpy as np
from math import prod
from nilearn.masking import (
        apply_mask,
        unmask
        )
from nilearn.plotting import (
        plot_epi,
        plot_stat_map,
        show)
from nilearn.image import (
        load_img,
        resample_to_img
        )
from nibabel import Nifti1Image
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
#import argparse

source_directory="/bcbl/home/public/MarcoMotion/Resting_State/analysis_timeSeries"
methods=['vanilla', 'nordic', 'tmmpca', 'mppca', 'nordic', 'nordic_magn']
subjects=['sub-001', 'sub-002', 'sub-003', 'sub-004', 'sub-005']
tasks=['task-HABLA1200', 'task-HABLA1700']
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

#parser=argparse.ArgumentParser(description="""Computes reliability for fMRI data over GM
#        so far this codes  get original data and split it in half""")
#parser.add_argument("--source_directory", default=None, type=str,
#        help="Full path to the source directory")
#parser.add_argument("--subjects", default=None, nargs="+",
#        help=""" subjects to iterate and do within method comparison
#        i.e. subjects=(sub-001, sub-002, sub-003)""")
#parser.add_argument("--tasks", default=None, nargs="+",
#        help=""" task to iterate and do within method comparison per task/run
#        i.e. tasks=(mppca,nordic,hydra,tmppca)""")
#parser.add_argument("--methods", default=None, nargs="+",
#        help=""" method to iterate and compare i.e. methods=(mppca,nordic,hydra,tmppca)""")
#parser.add_argument("--overwrite", default=True, type=bool,
#        help=""" Haults program if scatter plots exist. Default behaviour is True""")
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

#args = parser.parse_args()
#source_directory = args.source_directory
#subjects = args.subjects
#tasks = args.tasks
#methods = args.methods
#source_directory = args.source_directory
#overwrite = args.overwrite
#####################################################################################
###### Functions ####################################################################
#####################################################################################
#######################################################################################
#######################################################################################
def reliability_analysis(subject,task,methodx,mask,sbref,directory=source_directory,
        split=True,plot=False,savecorr=True,save=True,hist=True,residuals=False):

    print(f"""Worth double checking! To understand output """)
    print(f"""Computing reliability between halves of same process: {split} """)
    print(f"""Original time series is residual: {residuals}""")
    print(f"""Saving ... correlation matrixes: {savecorr}, r-values histogram: {hist}, plot: {save}""")
    ##############################
    print("Loading timeseries")
    ##############################
    part="OC"
    if residuals == True:
        part="residuals"
    
    if split == True:
        file1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+".nii.gz"
        print(f"""Loading epi file: {file1} while applying mask: {mask}""")
        epi_mask=apply_mask(file1,mask)
        shape=epi_mask.shape
        print(" Data loaded and masked!")
        print(f"""Mask: {mask} contains {shape[0]} voxels""")
        epi_mask=np.transpose(epi_mask)
        # splitting volume in two 
        epi_half1=epi_mask[:,:(shape[-1]//2)]
        epi_half2=epi_mask[:,(shape[-1]//2):]
        print(" Original epi time series divided in two")
    elif split == False:
        file1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+"1.nii.gz"
        file2=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+"2.nii.gz"
        print(f"""Loading epi file: {file1}""")
        epi_half1=apply_mask(file1,mask)
        shape1=epi_half1.shape
        print(f"""Loading epi file: {file2}""")
        epi_half2=apply_mask(file2,mask)
        shape2=epi_half2.shape
        if len(shape1) != len(shape2):
            raise ValueError(f"--ERROR-- data 1 has {len(shape1)} dimentions, while data 2 has {len(shape2)}")
        print(f"""Mask: {mask} contains {shape1[0]} voxels""")
        epi_half1=np.transpose(epi_half1)
        epi_half2=np.transpose(epi_half2)
    correlation_matrix_half1=np.corrcoef(epi_half1)
    print(f""" Functional connectivity for first half computed (pearson correlation) with shape {correlation_matrix_half1.shape}""")
    correlation_matrix_half2=np.corrcoef(epi_half2)
    print(f""" Functional connectivity for second half computed (pearson correlation) with shape {correlation_matrix_half2.shape}""")
    if savecorr == True:
        file_corr1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+split+"_half1.csv"
        np.savetxt(file_corr1,correlation_matrix_half1,delimiter=",")
        print(f""" Functional connectivity for first half saved as {file_corr1}""")
        file_corr2=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+split+"_half2.csv"
        np.savetxt(file_corr2,correlation_matrix_half2,delimiter=",")
        print(f""" Functional connectivity for second half saved as {file_corr2}""")
    #del epi_half1
    #del epi_half2
    #del epi_half1_mean
    #del epi_half2_mean
    reliability_vector=pearsonr(correlation_matrix_half1, correlation_matrix_half2).statistic
    if hist==True:
        fig, ax =plt.subplots(nrows=1,ncols=1)
        ax.hist(reliability_vector,bins=100,density=True,edgecolor='black')
        plt.xlabel("Coefficient values")
        plt.ylabel("Frequency")
        fig.suptitle("Reliability coefficients histogram")
        fig.savefig(directory+"/"+subject+"_ses-1_"+task+"_"+part+"_reliability_coefficients_"+methodx+split+"_histogram.png")
        plt.close(fig)
    print(f""" Reliability score computed voxel wise between two halves. result has shape: {reliability_vector.shape}""")
    if save == True:
        file_reliability=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+split+"_reliability.csv"
        np.savetxt(file_reliability,reliability_vector,delimiter=",")
        print(f""" Reliability vector saved as {file_reliability}""")
    if plot == True:
        plot_results=unmask(reliability_vector,mask)
        shape=plot_results.shape
        sbref_epi=load_img(sbref)
        #sbref_affine=sbref_epi.affine
        plot_results_affined=Nifti1Image(plot_results.get_fdata(),affine=sbref_epi.affine, header=sbref_epi.header)
        #plot_results_resampled=resample_to_img(plot_results,sbref_epi, interpolation='linear')
        #plot_results_resampled2=resample_img(plot_results,target_affine=sbref_epi.affine,target_shape=sbref_epi.shape,interpolation='nearest', force_resample=True,copy_header=False)
        print("Created new nilearn object to visualize results")
        #brain_reliability=plot_epi(plot_results_affined,bg_img=sbref_epi,colorbar=True,draw_cross=False,cut_coords=((shape[0]//2),(shape[1]//2),(shape[2]//2)),cmap="inferno",vmin=0,vmax=0.5)
        title=("Reliability map for "+subject+" "+methodx)
        brain_reliability=plot_stat_map(plot_results_affined,sbref_epi,colorbar=True,draw_cross=False,title=title,cut_coords=((shape[0]//2),(shape[1]//2),(shape[2]//2)),cmap="inferno",vmin=0,vmax=0.5)
        brain_reliability.savefig(directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+split+"_reliability.png")

#############################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            mask=source_directory+"/"+subject+"_ses-1_"+task+"_echo-1_part-mag_gm_mask-union.nii.gz"
            sbref="/bcbl/home/public/MarcoMotion/Resting_State/analysis/"+subject+"_ses-1_"+task+"_echo-1_part-mag_masked_sbref.nii.gz"
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,method,mask,sbref,plot=True)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,method,mask,sbref,plot=True, residuals=True)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,method,mask,sbref,split=False,plot=True)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
