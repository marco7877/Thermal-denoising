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
def reliability_analysis(subject,task,methodx,mask,directory=source_directory,
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
        file1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+"_gm_union.nii.gz"
        file2=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_gm_union.nii.gz"
        #file2=file2.replace("analysis_timeSeries","analysis")
        print(f"""Loading epi file: {file1}""")
        data_episeries_x=load_img(file1)
        data_episeries_x=np.array(data_episeries_x.dataobj)
        shape=data_episeries_x.shape
        data_mask=load_img(mask)
        data_mask=np.array(data_mask.dataobj)
        # creatind a 2D array of voxels * time series
        data_episeries_x=np.reshape(data_episeries_x,((prod(shape[0:-1])),shape[-1]))
        data_mask=np.reshape(data_mask,((prod(shape)),0))
        print("data reshaped")
        # creating a mask of non zero values so only getting GM 
        epi_mask=np.ma.masked_where(data_mask[:,0]!=0)
        print(f"""masked created for {sum(epi_mask.mask)} voxels""")
        # splitting volume in two 
        epi_half1=data_episeries_x[epi_mask.mask,:(shape[-1]//2)]
        epi_half2=data_episeries_x[epi_mask.mask,(shape[-1]//2):]
        print(" Original epi time series divided in two")
    elif split == False:
        file1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+"1_gm_union.nii.gz"
        file2=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_part-mag_bold_"+methodx+"2_gm_union.nii.gz"
        file3=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_gm_union.nii.gz"
        print(f"""Loading epi file: {file1}""")
        data_episeries_x1=load_img(file1)
        data_episeries_x1=np.array(data_episeries_x1.dataobj)
        shape1=data_episeries_x1.shape
        print(f"""Loading epi file: {file2}""")
        data_episeries_x2=load_img(file2)
        data_episeries_x2=np.array(data_episeries_x2.dataobj)
        shape2=data_episeries_x2.shape
        if len(shape1) != len(shape2):
            raise ValueError(f"--ERROR-- data 1 has {len(shape1)} dimentions, while data 2 has {len(shape2)}")
        # creatind a 2D array of voxels * time series
        data_episeries_x1=np.reshape(data_episeries_x1,((prod(shape1[0:-1])),shape1[-1]))
        data_episeries_x2=np.reshape(data_episeries_x2,((prod(shape2[0:-1])),shape2[-1]))
        print("data reshaped")
        data_mask=load_img(mask)
        data_mask=np.array(data_mask.dataobj)
        # creating a mask of non zero values so only getting GM 
        epi_mask=np.ma.masked_where(data_mask[:,0]!=0)
        print(f"""masked created for {sum(epi_mask.mask)} voxels""")
        epi_half1=data_episeries_x1[epi_mask.mask,:]
        epi_half2=data_episeries_x2[epi_mask.mask,:]
        #TODO: finish false condition where reliability gets computed between independent halves
    correlation_matrix_half1=np.corrcoef(epi_half1)
    print(f""" Functional connectivity for first half computed (pearson correlation) with shape {correlation_matrix_half1.shape}""")
    correlation_matrix_half2=np.corrcoef(epi_half2)
    print(f""" Functional connectivity for second half computed (pearson correlation) with shape {correlation_matrix_half2.shape}""")
    if savecorr == True:
        file_corr1=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+"_half1.csv"
        np.savetxt(file_corr1,correlation_matrix_half1,delimiter=",")
        print(f""" Functional connectivity for first half saved as {file_corr1}""")
        file_corr2=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+"_half2.csv"
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
        fig.savefig(directory+"/"+subject+"_ses-1_"+task+"_"+part+"_reliability_coefficients_"+methodx+"_histogram.png")
        plt.close(fig)
    print(f""" Reliability score computed voxel wise between two halves. result has shape: {reliability_vector.shape}""")
    if save == True:
        file_reliability=directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+"_reliability.csv"
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
        brain_reliability.savefig(directory+"/"+subject+"_ses-1_"+task+"_"+part+"_functional_connectivity_"+methodx+"_reliability.png")

#############################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            mask=directory+"/"+subject+"_ses-1_"+task+"_echo-1_part-mag_gm-alligned_mask.nii.gz"
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,mask,method,plot=True)
            except Reliability split halves:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,mask,method,plot=True, residuals=True)
            except Reliability residuals:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                reliability_analysis(subject,task,mask,method,split=False,plot=True)
            except Reliability diferent halves:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
