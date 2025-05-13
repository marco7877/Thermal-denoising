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
        compute_epi_mask
        )
from nilearn.plotting import (
        plot_epi,
        show)
from nilearn.image import(
        load_img,
        new_img_like
        )
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
def reliability_analysis(subject,task,methodx,directory=source_directory,
        plot=False,savecorr=False,save=True,hist=False):
    ##############################
    print("Loading timeseries")
    ##############################
    file1=directory+"/"+subject+"_ses-1_"+task+"_residuals_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    file2=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    #file2=file2.replace("analysis_timeSeries","analysis")
    print(f"""Loading epi file: {file1}""")
    data_episeries_x=load_img(file1)
    data_episeries_x=np.array(data_episeries_x.dataobj)
    data_episeries_x2=load_img(file2)
    data_episeries_x2=np.array(data_episeries_x2.dataobj)
    shape=data_episeries_x.shape
    # creatind a 2D array of voxels * time series
    data_episeries_x=np.reshape(data_episeries_x,((prod(shape[0:-1])),shape[-1]))
    data_episeries_x2=np.reshape(data_episeries_x2,((prod(shape[0:-1])),shape[-1]))
    print("data reshaped")
    # creating a mask of non zero values so only getting GM 
    epi_mask=np.ma.masked_where(data_episeries_x2[:,0]!=0,data_episeries_x2[:,0])
    print(f"""masked created for {sum(epi_mask.mask)} voxels""")
    # splitting volume in two 
    epi_half1=data_episeries_x[epi_mask.mask,:(shape[-1]//2)]
    epi_half2=data_episeries_x[epi_mask.mask,(shape[-1]//2):]
    print(" Original epi time series divided in two")
    # TODO: normalize voxelwise maybe pct change respect to mean
    #epi_half1_mean=np.mean(epi_half1,axis=1)
    #epi_half2_mean=np.mean(epi_half2,axis=1)
    #scaled_epi_half1=(np.subtract(epi_half1,epi_half1_mean[:,np.newaxis]))/epi_half1_mean[:,np.newaxis]
    #scaled_epi_half2=(np.subtract(epi_half2,epi_half2_mean[:,np.newaxis]))/epi_half2_mean[:,np.newaxis]
    print(" Halves re-scaled to represent pct change with respect to mean time series value per voxel")
    correlation_matrix_half1=np.corrcoef(epi_half1)
    print(f""" Functional connectivity for first half computed (pearson correlation) with shape {correlation_matrix_half1.shape}""")
    correlation_matrix_half2=np.corrcoef(epi_half2)
    print(f""" Functional connectivity for second half computed (pearson correlation) with shape {correlation_matrix_half2.shape}""")
    if savecorr == True:
        file_corr1=directory+"/"+subject+"_ses-1_"+task+"_residuals_functional_connectivity_"+methodx+"_half1.csv"
        np.savetxt(file_corr1,correlation_matrix_half1,delimiter=",")
        print(f""" Functional connectivity for first half saved as {file_corr1}""")
        file_corr2=directory+"/"+subject+"_ses-1_"+task+"_residuals_functional_connectivity_"+methodx+"_half2.csv"
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
        fig.savefig(directory+"/"+subject+"_ses-1_"+task+"_residuals_reliability_coefficients_"+methodx+"_histogram.png")
        plt.close(fig)
    print(f""" Reliability score computed voxel wise between two halves. result has shape: {reliability_vector.shape}""")
    if save == True:
        file_reliability=directory+"/"+subject+"_ses-1_"+task+"_residuals_functional_connectivity_"+methodx+"_reliability.csv"
        np.savetxt(file_reliability,reliability_vector,delimiter=",")
        print(f""" Reliability vector saved as {file_reliability}""")
    reliability_target=np.zeros((prod(shape[:-1])))
    reliability_target[epi_mask.mask]=reliability_vector**2
    del reliability_vector
    if plot == True:
        brain_space=np.reshape(reliability_target,(shape[:-1]))
        print(f""" Moved results back to original space""")
        ref_img=load_img(file2)
        plot_results=new_img_like(ref_img,brain_space)
        print("Created new nilearn object to visualize results")
        brain_reliability=plot_epi(plot_results,colorbar=True,draw_cross=False,cut_coords=((shape[0]//2),(shape[1]//2),(shape[2]//2)),cmap="inferno",vmin=0,vmax=0.5)
        brain_reliability.savefig(directory+"/"+subject+"_ses-1_"+task+"_residuals_functional_connectivity_"+methodx+"_reliability.png")

#############################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,method,plot=True,hist=True)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
