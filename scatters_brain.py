#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

from nilearn import image
from nilearn import signal
from nilearn.masking import (
        apply_mask,
        compute_epi_mask
        )
from scipy import stats
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
from itertools import combinations
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
    epi_file1=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    epi_file2=directory+"/"+subject+"_ses-1_"+task+"_OC_part-mag_bold_"+methody+"_gm_union.nii.gz"
    if automask == True:
        # We assume both epi are on the same space and mask
        mask=compute_epi_mask(epi_file1,exclude_zeros=True)
    else:
        mask=compute_epi_mask(mask_file,exclude_zeros=True)
    masked_datay=apply_mask(epi_file1,mask)
    masked_datax=apply_mask(epi_file1,mask)
    df_plot=pd.DataFrame({'X':masked_datax,'Y':masked_datay})
    plotnine=(ggplot(df_plot,aes("X","Y"))
            +geom_point(alpha=0.3,colour="gray")
            +geom_smooth(method="lm")
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,150])
            +scale_y_continuous(limits=[0,150])
            +xlab(task+" "+methodx)
            +ylab(task+" "+methody)
            +theme(
                panel_background=element_rect(fill="white"),
                panel_grid=element_blank(),
                panel_border=element_blank())
            )
    plotnine.draw()
    plotnine.save(directory+"/"+subject+task+"vanilla-"+methody+".png",verbose=False)
#####################################################################################
###### Main      ####################################################################
#####################################################################################

