#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

from nilearn import image
from nilearn import signal
from nilearn.masking import apply_mask
from nilearn.maskers import NiftiMasker
from scipy import stats
import numpy as np
import pandas as pd
import os
import argparse
from plotnine import  *

#####################################################################################
###### Arguments ####################################################################
#####################################################################################

parser=argparse.ArgumentParser(description="""Generates scatterplots with plotnine
        from fMRI data with only one brik """)
parser.add_argument("--source_directory", default=None, type=str,
        help="Full path to the BIDS directory")
parser.add_argument("--methods", default=None, type=list,
        help=""" method to iterate and compare i.e. methods=(mppca,nordic,hydra,tmppca)""")
parser.add_argument("--nii_pattern", default=None, type=str,
        help="""file pattern to search with subject, task and method replaced with
        SUB, TSK, and MTHD. i.e. SUB_ses-1_TSK_OC_part-mag_bold_MTHD_gm_union""")
#####################################################################################
###### Arguments ####################################################################
#####################################################################################

args = parser.parse_args()
bids_dir = args.bids_dir
mprageize_dir = args.mprageize_dir
mprageize = mprageize_dir + "/3dMPRAGEize"
overwrite = args.overwrite

#parser.add_argument("--tsnr_extention", default=None, type=str,
#                    help="""Extention of the target tsnr file to be look upon
#                    i.e., if my target file has the name:
#                    sub-001_task-HABLA1200_masked_epi_gm_ocDenoised_tsnr.nii.gz
#                    then, OC_extention=task-HABLA1200_masked_epi_gm_ocDenoised_tsnr.nii.gz""")
#
#source_directory = "/bcbl/home/public/MarcoMotion/Resting_State/analysis/"
#task="HABLA1200"
#tsnr_class="1200-OHBM"
#tsnr_extention="task-"+task+"_OC_part-mag_bold_"

#ext="_gm.nii.gz"
#methods=["hydra", "nordic","vanilla", "tmmpca", "mppca", "nordic_magn"]
#labels=["hydra", "nordic","vanilla", "tmmpca", "mppca", "nordic_magn"]
#voxel="2.4*2.4*2.4mm"
