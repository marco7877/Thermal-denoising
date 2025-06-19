#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""
import numpy as np
#from math import prod
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
from itertools import combinations
#import os
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
#def reliability_analysis(subject,task,methodx,mask,sbref,directory=source_directory,
def reliability_analysis(files,mask,sbref,directory=source_directory,
        split=True,plot=False,savecorr=False,hist=True,make_nifti=False):

    print(f"""Worth double checking! To understand output """)
    print(f"""Computing reliability between halves of same process: {split} """)
    print(f"""Saving ... correlation matrixes: {savecorr}, r-values histogram: {hist}, plot: {plot}""")
    ##############################
    print("Loading timeseries")
    ##############################
    array_dict={}
    for i in list(range(len(files))):
        print(f"""Loading epi file: {files[i]} while applying mask: {mask}""")
        array_dict[i]=np.transpose(apply_mask(files[i],mask))
        print(" Data loaded and masked!")
        shape=array_dict[i].shape
        print(f"""Mask: {mask} contains {shape[0]} voxels""")
    if len(array_dict) == 1:
        corr_dict={}
        corr_dict[0]=np.corrcoef(array_dict[0][:,:(shape[-1]//2)])
        print(f""" Functional connectivity for computed (pearson correlation) with shape {corr_dict[0].shape}""")
        corr_dict[1]=np.corrcoef(array_dict[0][:,(shape[-1]//2):])
        print(f""" Functional connectivity for computed (pearson correlation) with shape {corr_dict[1].shape}""")
        print(" Original epi time series divided in two")
        files.append(files[0])
        for i in range(2):
            files[i].replace(files[i].split("_")[-1].split(".")[0],files[i].split("_")[-1].split(".")[0]+str(i))
            if savecorr == True:
                np.savetxt(files[i].replace(files[i].split("_")[-1],files[i].split("_")[-1].split(".")[0]+"fconnectivity.csv"),corr_dict[i],delimiter=",")
                print(f""" Functional connectivity for saved as {files[i].replace(files[i].split("_")[-1],files[i].split("_")[-1].split(".")[0]+"fconnectivity.csv")}""")
    elif len(array_dict) > 1:
        for i in list(range(len(files))):
            corr_dict[i]=np.corrcoef(array_dict[i])
            print(f""" Functional connectivity for computed (pearson correlation) with shape {corr_dict[i].shape}""")
            if savecorr == True:
                np.savetxt(files[i].replace(files[i].split("_")[-1],files[i].split("_")[-1].split(".")[0]+"fconnectivity.csv"),corr_dict[i],delimiter=",")
                print(f""" Functional connectivity for saved as {files[i].replace(files[i].split("_")[-1],files[i].split("_")[-1].split(".")[0]+"fconnectivity.csv")}""")
    perm_volumes=list(combinations(list(range(len(corr_dict))),2))
    print(f""" Calculating reliability for combinations""")
    reliability_dict={}
    for i in range(len(perm_volumes)):
        reliability_dict[i]=pow(pearsonr(corr_dict[perm_volumes[i][0]],corr_dict[perm_volumes[i][1]]).statistic,2)
        print(f"""Reliability calculated for epi combinaiton {1+i}""")
        if make_nifti == True:
            plot_results=unmask(reliability_dict[i],mask)
            plot_results.to_filename(files[perm_volumes[i][0]].replace(files[perm_volumes[i][0]].split("_")[-1].split(".")[0],files[perm_volumes[i][0]].split("_")[-1].split(".")[0]+str(perm_volumes[i][1])+".nii.gz"))
        if hist == True:
            fig, ax =plt.subplots(nrows=1,ncols=1)
            ax.hist(reliability_dict[i],bins=100,density=True,edgecolor='black')
            plt.xlabel("Coefficient values")
            plt.ylabel("Frequency")
            fig.suptitle("Reliability coefficients histogram")
            fig.savefig(files[i].replace(files[perm_volumes[i][0]].split("_")[-1],files[perm_volumes[i][0]].split("_")[-1].split(".")[0]+"_histogram.png"))
            plt.close(fig)

        if plot == True:
            plot_results=unmask(reliability_dict[i],mask)
            shape_epi=plot_results.shape
            sbref_epi=load_img(sbref)
            print("Loaded sbref for background: {sbref}")
            plot_results_affined=Nifti1Image(plot_results.get_fdata(),affine=sbref_epi.affine, header=sbref_epi.header)
            print("Created new nilearn object to visualize results")
            title=("Reliability map for "+files[perm_volumes[i][0]].split("_")[0].split("/")[-1]+" "+files[perm_volumes[i][0]].split("_")[-1].split(".")[0])
            brain_reliability=plot_stat_map(plot_results_affined,sbref_epi,colorbar=True,draw_cross=False,title=title,cut_coords=((shape_epi[0]//2),(shape_epi[1]//2),(shape_epi[2]//2)),cmap="inferno",vmin=0,vmax=0.5)
            brain_reliability.savefig(files[perm_volumes[i][0]].replace(files[perm_volumes[i][0]].split("_")[-1],files[perm_volumes[i][0]].split("_")[-1].split(".")[0]+"_reliability.png"))

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
                reliability_analysis(subject,task,method,mask,sbref,hist=False,plot=True)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
            try:
                print(f"""############################################################################""")
                print(f"""##########reliability_analysis({subject},{task},{method},{mask},{sbref},plot=True,residuals=True)
#######################""")
                reliability_analysis(subject,task,method,mask,sbref,plot=True,residuals=True)
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
