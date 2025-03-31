#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

from nilearn.masking import (
        apply_mask,
        compute_epi_mask
        )
import numpy as np
import pandas as pd
import os
import argparse
from plotnine import *
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

parser=argparse.ArgumentParser(description="""Generates scatterplots with plotnine
        from fMRI data with only one brik """)
parser.add_argument("--source_directory", default=None, type=str,
        help="Full path to the source directory")
parser.add_argument("--subjects", default=None, nargs="+",
        help=""" subjects to iterate and do within method comparison
        i.e. subjects=(sub-001, sub-002, sub-003)""")
parser.add_argument("--tasks", default=None, nargs="+",
        help=""" task to iterate and do within method comparison per task/run
        i.e. tasks=(mppca,nordic,hydra,tmppca)""")
parser.add_argument("--methods", default=None, nargs="+",
        help=""" method to iterate and compare i.e. methods=(mppca,nordic,hydra,tmppca)""")
parser.add_argument("--overwrite", default=True, type=bool,
        help=""" Haults program if scatter plots exist. Default behaviour is True""")
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

args = parser.parse_args()
source_directory = args.source_directory
subjects = args.subjects
tasks = args.tasks
methods = args.methods
source_directory = args.source_directory
overwrite = args.overwrite
#####################################################################################
###### Functions ####################################################################
#####################################################################################
def scatter_plot(subject,task,methodx,methody,directory=source_directory,automask=True,mask_file=""):
    gm_file1=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    gm_file2=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methody+"_gm_union.nii.gz"
    ngm_file1=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_non-gm_union.nii.gz"
    ngm_file2=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methody+"_non-gm_union.nii.gz"
    if automask == True:
        # We assume both epi are on the same space and mask
        gm_mask=compute_epi_mask(gm_file1,exclude_zeros=True)
        ngm_mask=compute_epi_mask(ngm_file1,exclude_zeros=True)
    else:
        gm_mask=compute_epi_mask(mask_file,exclude_zeros=True)
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=apply_mask(gm_file2,gm_mask)
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=apply_mask(gm_file1,gm_mask)
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=apply_mask(ngm_file2,ngm_mask)
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=apply_mask(ngm_file1,ngm_mask)
    df_gm=pd.DataFrame({'X':masked_data_gm_x,'Y':masked_data_gm_y,"type":"gray matter"})
    df_ngm=pd.DataFrame({'X':masked_data_ngm_x,'Y':masked_data_ngm_y,"type":"non-gray matter"})
    df_plot=pd.concat([df_gm,df_ngm],axis=0)
    df_plot=df_plot.replace(0,np.nan)
    df_plot=df_plot.dropna(axis=0)
    plotnine_img=(ggplot(df_plot,aes("X","Y"))
            +geom_point(alpha=0.6)
#            +geom_smooth(color=["black","black"])
#            +geom_smooth(method="lm")
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,200])
            +scale_y_continuous(limits=[0,400])
            +xlab("tSNR "+methodx)
            +ylab("tSNR "+methody)
            +theme_classic())
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+".png",verbose=False)
#######################################################################################
def scatter_plotS0T2(subject,task,methodx,methody,directory=source_directory,automask=True,mask_file=""):
    directory=directory.replace("analysis","T2")
    ##############################
    print("T2")
    ##############################
    gm_file1=directory+"/"+subject+"_ses-1_"+task+"_T2_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    gm_file2=directory+"/"+subject+"_ses-1_"+task+"_T2_part-mag_bold_"+methody+"_gm_union.nii.gz"
    ngm_file1=directory+"/"+subject+"_ses-1_"+task+"_T2_part-mag_bold_"+methodx+"_non-gm_union.nii.gz"
    ngm_file2=directory+"/"+subject+"_ses-1_"+task+"_T2_part-mag_bold_"+methody+"_non-gm_union.nii.gz"
    if automask == True:
        # We assume both epi are on the same space and mask
        gm_mask=compute_epi_mask(gm_file1,exclude_zeros=True)
        ngm_mask=compute_epi_mask(ngm_file1,exclude_zeros=True)
    else:
        gm_mask=compute_epi_mask(mask_file,exclude_zeros=True)
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=apply_mask(gm_file2,gm_mask)
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=apply_mask(gm_file1,gm_mask)
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=apply_mask(ngm_file2,ngm_mask)
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=apply_mask(ngm_file1,ngm_mask)
    df_gm=pd.DataFrame({'X':masked_data_gm_x,'Y':masked_data_gm_y,"type":"gray matter"})
    df_ngm=pd.DataFrame({'X':masked_data_ngm_x,'Y':masked_data_ngm_y,"type":"non-gray matter"})
    df_plot=pd.concat([df_gm,df_ngm],axis=0)
    df_plot.loc[df_plot["X"]>0.2,"X"]=np.nan
    df_plot=df_plot.replace(0,np.nan)
    df_plot=df_plot.dropna(axis=0)
    plotnine_img=(ggplot(df_plot,aes("X","Y"))
            +geom_point(alpha=0.9)
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,0.2])
            +scale_y_continuous(limits=[0,0.2])
            +xlab("T2* "+methodx)
            +ylab("T2* "+methody)
            +theme_classic())
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+"_T2.png",verbose=False)
    ##############################
    print("S0")
    ##############################
    gm_file1=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    gm_file2=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methody+"_gm_union.nii.gz"
    ngm_file1=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methodx+"_non-gm_union.nii.gz"
    ngm_file2=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methody+"_non-gm_union.nii.gz"
    if automask == True:
        # We assume both epi are on the same space and mask
        gm_mask=compute_epi_mask(gm_file1,exclude_zeros=True)
        ngm_mask=compute_epi_mask(ngm_file1,exclude_zeros=True)
    else:
        gm_mask=compute_epi_mask(mask_file,exclude_zeros=True)
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=apply_mask(gm_file2,gm_mask)
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=apply_mask(gm_file1,gm_mask)
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=apply_mask(ngm_file2,ngm_mask)
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=apply_mask(ngm_file1,ngm_mask)
    df_gm=pd.DataFrame({'X':masked_data_gm_x,'Y':masked_data_gm_y,"type":"gray matter"})
    df_ngm=pd.DataFrame({'X':masked_data_ngm_x,'Y':masked_data_ngm_y,"type":"non-gray matter"})
    df_plot=pd.concat([df_gm,df_ngm],axis=0)
    df_plot=df_plot.replace(0,np.nan)
    df_plot=df_plot.dropna(axis=0)
    plotnine_img=(ggplot(df_plot,aes("X","Y"))
            +geom_point(alpha=0.9)
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,200000])
            +scale_y_continuous(limits=[0,200000])
            +xlab("S0 "+methodx)
            +ylab("S0 "+methody)
            +theme_classic())   
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+"_S0.png",verbose=False)
#######################################################################################
def scatter_plotT2sPCT(subject,task,methodx,directory=source_directory,automask=True,mask_file=""):
    ##############################
    print("T2* series (%)")
    ##############################
    directory=directory.replace("analysis","T2")
    file1=directory+"/"+subject+"_ses-1_"+task+"_T2spc_part-mag_bold_"+methodx+".nii.gz"
    if automask == True:
        # We assume both epi are on the same space and mask
        mask=compute_epi_mask(file1,exclude_zeros=True)
    else:
        gm_mask=compute_epi_mask(mask_file,exclude_zeros=True)
    print(f"""Loading T2 series file: {file1}""")
    masked_data_t2s_x=apply_mask(file1,mask)
    shape=masked_data_t2s_x.shape
    masked_data_t2s_x=np.reshape(masked_data_t2s_x,(shape[0]*shape[1]))
    df_t2s=pd.DataFrame({'X':masked_data_t2s_x,"type":"T2*"})
    df_t2s.loc[df_t2s["X"]>0.01,"X"]=np.nan
    df_t2s.loc[df_t2s["X"]< -0.01,"X"]=np.nan
    df_t2s=df_t2s.replace(0,np.nan)
    df_t2s=df_t2s.dropna(axis=0)
    plotnine_img=(ggplot(df_t2s,aes(x="X",y=after_stat("count/np.sum(count)")))
            +geom_histogram()
            +scale_x_continuous(limits=[-0.011,0.011])
            +xlab("T2* change(%)"+methodx)
            +ylab("Count (normalized) "+methodx)
            +theme_classic())
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+"_T2%.png",verbose=False)
    ##############################
    print("T2* series (CDF)")
    ##############################
    count,bins_count=np.histogram(df_t2s["X"],bins=100)
    prob_densityfunction=count/sum(count)
    cumulative_densityfuction=np.cumsum(prob_densityfunction)
    df_cdf=pd.DataFrame({"X":bins_count[1:],"Y":cumulative_densityfuction})
    plotnine_img=(ggplot(df_cdf,aes(x="X",y="Y"))
            +geom_line(linetype="dotted")
            +scale_x_continuous(limits=[-0.011,0.011])
            +scale_y_continuous(limits=[0,1])
            +xlab("T2* "+methodx)
            +ylab("Probability (cdf)")
            +theme_classic()
            )
    plotnine_img.save(directory+"/"+subject+task+"self-"+methody+"_T2(cdf).png",verbose=False)
#############################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            try:
                scatter_plot(subject,task,"vanilla",method)
                scatter_plotS0T2(subject,task,"vanilla",method)
                scatter_plotT2sPCT(subject,task,method)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
