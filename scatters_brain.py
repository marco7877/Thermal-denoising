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
from plotnine import (
        ggplot,
        aes,
        geom_point,
        scale_x_continuous,
        scale_y_continuous,
        xlab,
        ylab,
        geom_smooth,
        geom_abline,
        theme,
        panel_background,
        panel_grid,
        panel_border,
        element_blank,
        element_rect)
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

parser=argparse.ArgumentParser(description="""Generates scatterplots with plotnine
        from fMRI data with only one brik """)
parser.add_argument("--source_directory", default=None, type=str,
        help="Full path to the source directory")
parser.add_argument("--subjects", default=None, type=list,
        help=""" subjects to iterate and do within method comparison
        i.e. subjects=(sub-001, sub-002, sub-003)""")
parser.add_argument("--tasks", default=None, type=list,
        help=""" task to iterate and do within method comparison per task/run
        i.e. tasks=(mppca,nordic,hydra,tmppca)""")
parser.add_argument("--methods", default=None, type=list,
        help=""" method to iterate and compare i.e. methods=(mppca,nordic,hydra,tmppca)""")
parser.add_argument("--nii_pattern", default=None, type=list,
        help="""list of ordered strings with the elements to iterate over. They 
        must follow the pattern ([SUB]_ses-1,[_TSK]_OC_part-mag_bold,[_MTHD]_gm_union)
        where whats within square braces are the ommited variables we will be iterating""")
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
nii_pattern= args.nii_pattern
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
    plotnine_img=(ggplot(df_plot,aes("X","Y",color="type"))
            +geom_point(alpha=0.3)
            +geom_smooth(method="lm")
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,200])
            +scale_y_continuous(limits=[0,400])
            +xlab("tSNR "+methodx)
            +ylab("tSNR "+methody)
            +theme_classic())
    plotnine_img.draw()
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+".png",verbose=False)
#####################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            try:
                scatter_plot(subject,task,"vanilla",method)
            except:
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
