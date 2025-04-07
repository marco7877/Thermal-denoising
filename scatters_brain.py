#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""
from math import prod
from nilearn.masking import (
        apply_mask,
        compute_epi_mask
        )
from nilearn.image import(
        load_img
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
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=load_img(gm_file2)
    masked_data_gm_y=np.array(masked_data_gm_y.dataobj)
    masked_data_gm_y=np.reshape(masked_data_gm_y,(prod(masked_data_gm_y.shape)))
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=load_img(gm_file1)
    masked_data_gm_x=np.array(masked_data_gm_y.dataobj)
    masked_data_gm_x=np.reshape(masked_data_gm_x,(prod(masked_data_gm_x.shape)))
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=load_img(ngm_file2)
    masked_data_ngm_y=np.array(masked_data_ngm_y.dataobj)
    masked_data_ngm_y=np.reshape(masked_data_ngm_y,(prod(masked_data_ngm_y.shape)))
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=load_img(ngm_file1)
    masked_data_ngm_x=np.array(masked_data_ngm_x.dataobj)
    masked_data_ngm_x=np.reshape(masked_data_ngm_x,(prod(masked_data_ngm_x.shape)))
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
def scatter_plotS0R2(subject,task,methodx,methody,directory=source_directory):
    directory=directory.replace("analysis","T2")
    ##############################
    print("R2")
    ##############################
    gm_file1=directory+"/"+subject+"_ses-1_"+task+"_R2_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    gm_file2=directory+"/"+subject+"_ses-1_"+task+"_R2_part-mag_bold_"+methody+"_gm_union.nii.gz"
    ngm_file1=directory+"/"+subject+"_ses-1_"+task+"_R2_part-mag_bold_"+methodx+"_non-gm_union.nii.gz"
    ngm_file2=directory+"/"+subject+"_ses-1_"+task+"_R2_part-mag_bold_"+methody+"_non-gm_union.nii.gz"
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=load_img(gm_file2)
    masked_data_gm_y=np.array(masked_data_gm_y.dataobj)
    masked_data_gm_y=np.reshape(masked_data_gm_y,(prod(masked_data_gm_y.shape)))
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=load_img(gm_file1)
    masked_data_gm_x=np.array(masked_data_gm_x.dataobj)
    masked_data_gm_x=np.reshape(masked_data_gm_x,(prod(masked_data_gm_x.shape)))
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=load_img(ngm_file2)
    masked_data_ngm_y=np.array(masked_data_ngm_y.dataobj)
    masked_data_ngm_y=np.reshape(masked_data_ngm_y,(prod(masked_data_ngm_y.shape)))
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=load_img(ngm_file1)
    masked_data_ngm_x=np.array(masked_data_ngm_x.dataobj)
    masked_data_ngm_x=np.reshape(masked_data_ngm_x,(prod(masked_data_ngm_x.shape)))
    df_gm=pd.DataFrame({'X':masked_data_gm_x,'Y':masked_data_gm_y,"type":"gray matter"})
    df_ngm=pd.DataFrame({'X':masked_data_ngm_x,'Y':masked_data_ngm_y,"type":"non-gray matter"})
    df_plot=pd.concat([df_gm,df_ngm],axis=0)
    df_plot.loc[df_plot["X"]>50,"X"]=np.nan
    df_plot=df_plot.replace(0,np.nan)
    df_plot=df_plot.dropna(axis=0)
    plotnine_img=(ggplot(df_plot,aes("X","Y"))
            +geom_point(alpha=0.9)
            +geom_abline(intercept=0,slope=1,colour="red",linetype="dotted")
            +scale_x_continuous(limits=[0,50])
            +scale_y_continuous(limits=[0,50])
            +xlab("R2* "+methodx)
            +ylab("R2* "+methody)
            +theme_classic())
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methody+"_R2.png",verbose=False)
    ##############################
    print("S0")
    ##############################
    gm_file1=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methodx+"_gm_union.nii.gz"
    gm_file2=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methody+"_gm_union.nii.gz"
    ngm_file1=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methodx+"_non-gm_union.nii.gz"
    ngm_file2=directory+"/"+subject+"_ses-1_"+task+"_S0_part-mag_bold_"+methody+"_non-gm_union.nii.gz"
    print(f"""Loading gm file: {gm_file2} as y axis""")
    masked_data_gm_y=load_img(gm_file2)
    masked_data_gm_y=np.array(masked_data_gm_y.dataobj)
    masked_data_gm_y=np.reshape(masked_data_gm_x,(prod(masked_data_gm_y.shape)))
    print(f"""Loading gm file: {gm_file1} as x axis""")
    masked_data_gm_x=load_img(gm_file1,gm_mask)
    masked_data_gm_x=np.array(masked_data_gm_x.dataobj)
    masked_data_gm_x=np.reshape(masked_data_gm_x,(prod(masked_data_gm_x.shape)))
    print(f"""Loading non-gm file: {ngm_file2} as y axis""")
    masked_data_ngm_y=load_img(ngm_file2,ngm_mask)
    masked_data_ngm_y=np.array(masked_data_ngm_y.dataobj)
    masked_data_ngm_y=np.reshape(masked_data_ngm_y,(prod(masked_data_ngm_y.shape)))
    print(f"""Loading gm file: {ngm_file1} as x axis""")
    masked_data_ngm_x=load_img(ngm_file1,ngm_mask)
    masked_data_ngm_x=np.array(masked_data_ngm_x.dataobj)
    masked_data_ngm_x=np.reshape(masked_data_ngm_x,(prod(masked_data_ngm_x.shape)))
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
def scatter_plotR2sPCT(subject,task,methodx,directory=source_directory,automask=True,mask_file=""):
    ##############################
    print("R2* series (%)")
    ##############################
    directory=directory.replace("analysis","T2")
    file1=directory+"/"+subject+"_ses-1_"+task+"_R2spc_part-mag_bold_"+methodx+".nii.gz"
    print(f"""Loading R2 series file: {file1}""")
    masked_data_t2s_x=load_img(file1)
    masked_data_t2s_x=np.array(masked_data_t2s_x.dataobj)
    shape=masked_data_t2s_x.shape
    masked_data_t2s_x=np.reshape(masked_data_t2s_x,((prod(shape))))
    print("data reshaped")
    df_t2s=pd.DataFrame({'X':masked_data_t2s_x,"type":"R2*"})
    df_t2s.loc[df_t2s["X"]>=1,"X"]=np.nan
    df_t2s.loc[df_t2s["X"]<= -1,"X"]=np.nan
    df_t2s=df_t2s.replace(0,np.nan)
    df_t2s=df_t2s.dropna(axis=0)
    print("creating plotnine")
    plotnine_img=(ggplot(df_t2s,aes(x="X",y=after_stat("count/np.sum(count)")))
            +geom_histogram(fill="lightgray",alpha=0.7)
            +scale_x_continuous(limits=[-1.1,1.1])
            +xlab("Delta R2*"+methodx)
            +ylab("Count (normalized) "+methodx)
            +theme_classic())
    print("plotnine object ready")
    plotnine_img.save(directory+"/"+subject+task+"vanilla-"+methodx+"_R2%.png")
    ##############################
    print("R2* series (CDF)")
    ##############################
    count,bins_count=np.histogram(df_t2s["X"],bins=100)
    print("Calculating PDF")
    prob_densityfunction=count/sum(count)
    print("Calculating CDF")
    cumulative_densityfuction=np.cumsum(prob_densityfunction)
    df_cdf=pd.DataFrame({"X":bins_count[1:],"Y":cumulative_densityfuction})
    print("creating plotnine object")
    plotnine_img=(ggplot(df_cdf,aes(x="X",y="Y"))
            +geom_line(color="red",linetype="solid",size=1)
            +scale_x_continuous(limits=[-1.1,1.1])
            +scale_y_continuous(limits=[0,1])
            +xlab("R2* "+methodx)
            +ylab("Probability (cdf)")
            +theme_classic()
            )
    print("plotnine object created")
    plotnine_img.save(directory+"/"+subject+task+"self-"+methodx+"_R2(cdf).png",verbose=False)
#############################################################################################
###### Main      ####################################################################
#####################################################################################
for subject in subjects:
    for task in tasks:
        for method in methods:
            try:
                print(f"""############################################################################""")
                print(f"""##########scatter_plot({subject},{task},"vanilla",{method})#######################""")
                #scatter_plot(subject,task,"vanilla",method)
                print(f"""############################################################################""")
                print(f"""##########scatter_plotS0R2({subject},{task},"vanilla",{method})#######################""")
                scatter_plotS0R2(subject,task,"vanilla",method)
                print(f"""############################################################################""")
                print(f"""##########scatter_plotR2sPCT({subject},{task},{method})#######################""")
                #scatter_plotR2sPCT(subject,task,method)
            except:
                print(f"""############################################################################""")
                print(f"""############################################################################""")
                print(f"""########################  ERROR  #############  ERROR  #####################""")
                print(f"""Something went wrong for subject: {subject}, task:{task}, and method:{method}""")
                print(f"""############################################################################""")
                print(f"""############################################################################""")
