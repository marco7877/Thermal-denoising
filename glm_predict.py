#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023

@author: mflores
"""

from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from nibabel import Nifti1Image
from nilearn.image import load_img, resample_to_img
from nilearn.masking import apply_mask, unmask
from nilearn.plotting import plot_epi, plot_stat_map, show
from scipy.stats import pearsonr

# import os
# import argparse

###############################################
###### Arguments ##############################
###############################################

# parser=argparse.ArgumentParser(description="Computes reliability for fMRI data over GM
#        so far this codes  get original data and split it in half")
# parser.add_argument("--source_dir", default=None, type=str,
#        help="Full path to the source directory")
# parser.add_argument("--subjects", default=None, nargs="+",
#        help=" subjects to iterate and do within method comparison
#        i.e. subjects=(sub-001, sub-002, sub-003)")
# parser.add_argument("--tasks", default=None, nargs="+",
#        help=" task to iterate and do within method comparison per task/run
#        i.e. tasks=(mppca,nordic,hydra,tmppca)")
# parser.add_argument("--methods", default=None, nargs="+",
#        help=" method to iterate and compare i.e. methods=(mppca,nordic,hydra,tmppca)")
# parser.add_argument("--overwrite", default=True, type=bool,
#        help=" Haults program if scatter plots exist. Default behaviour is True")
###############################################
###### Arguments ##############################
###############################################

# args = parser.parse_args()
# source_dir = args.source_dir
# subjects = args.subjects
# tasks = args.tasks
# methods = args.methods
# source_dir = args.source_dir
# overwrite = args.overwrite
###############################################
###### Functions ##############################
###############################################
def reliability_analysis(
    epi_fname, mask, sbref, plot=True, savecorr=False, hist=True, make_nifti=True
):
    # Define variables
    array_dict = {}
    corr_dict = {}
    reliability_dict = {}

    print("Worth double checking! To understand output ")
    print(
        f"Saving ... correlation matrixes: {savecorr}, r-values histogram: {hist}, plot: {plot}"
    )
    ##############################
    print("Loading timeseries")
    ##############################
    for i in list(range(len(epi_fname))):
        print(f"Loading epi file: {epi_fname[i]} while applying mask: {mask}")
        array_dict[i] = np.transpose(apply_mask(epi_fname[i], mask))
        print(" Data loaded and masked!")
        shape = array_dict[i].shape
        print(f"Mask: {mask} contains {shape[0]} voxels")

    if len(array_dict) == 1:
        # Single run, split in two
        corr_dict[0] = np.corrcoef(array_dict[0][:, : (shape[-1] // 2)])
        print(
            f"Functional connectivity for computed (pearson correlation) with shape {corr_dict[0].shape}"
        )
        corr_dict[1] = np.corrcoef(array_dict[0][:, (shape[-1] // 2) :])
        print(
            f"Functional connectivity for computed (pearson correlation) with shape {corr_dict[1].shape}"
        )
        print(" Original epi time series divided in two")
        perm_volumes = list(combinations(range(len(corr_dict)), 2))
        epi_fname.append(epi_fname[0])

        for i in range(2):
            base_fname = epi_fname[i].split("_")[-1]
            # TODO what is this?
            epi_fname[i].replace(
                base_fname.split(".")[0],
                base_fname.split(".")[0] + str(i),
            )
            if savecorr:
                # TODO perm_volumes is not defined before we hit this point?
                np.savetxt(
                    epi_fname[i].replace(
                        base_fname,
                        base_fname.split(".")[0]
                        + str(perm_volumes[i][0])
                        + str(perm_volumes[i][1])
                        + "fconnectivity.csv",
                    ),
                    corr_dict[i],
                    delimiter=",",
                )
                print(
                    f"Functional connectivity for saved as {epi_fname[i].replace(epi_fname[i].split('_')[-1], epi_fname[i].split('_')[-1].split('.')[0] + 'fconnectivity.csv')}"
                )

    elif len(array_dict) > 1:
        # Multiple runs, compute correlation for each
        perm_volumes = list(combinations(range(len(array_dict)), 2))
        for i in list(range(len(epi_fname))):
            base_fname = epi_fname[i].split("_")[-1]
            corr_dict[i] = np.corrcoef(array_dict[i])
            print(
                f"Functional connectivity for computed (pearson correlation) with shape {corr_dict[i].shape}"
            )
            # TODO should this be corr_dict[i]
            if savecorr:
                np.savetxt(
                    epi_fname[i].replace(
                        base_fname,
                        base_fname.split(".")[0]
                        + str(perm_volumes[i][0])
                        + str(perm_volumes[i][1])
                        + "fconnectivity.csv",
                    ),
                    corr_dict[i],
                    delimiter=",",
                )
                print(
                    f"Functional connectivity for saved as {epi_fname[i].replace(epi_fname[i].split('_')[-1], epi_fname[i].split('_')[-1].split('.')[0] + 'fconnectivity.csv')}"
                )

    print(" Calculating reliability for combinations")
    for i in range(len(perm_volumes)):
        reliability_dict[i] = pow(
            pearsonr(
                corr_dict[perm_volumes[i][0]], corr_dict[perm_volumes[i][1]]
            ).statistic,
            2,
        )
        print(f"Reliability calculated for epi combinaiton {1 + i}")
        if make_nifti:
            plot_results = unmask(reliability_dict[i], mask)
            plot_results.to_filename(
                epi_fname[perm_volumes[i][0]].replace(
                    epi_fname[perm_volumes[i][0]].split("_")[-1].split(".")[0],
                    epi_fname[perm_volumes[i][0]].split("_")[-1].split(".")[0]
                    + str(perm_volumes[i][0])
                    + str(perm_volumes[i][1]),
                )
            )

        if hist:
            fig, ax = plt.subplots(nrows=1, ncols=1)
            ax.hist(reliability_dict[i], bins=100, density=True, edgecolor="black")
            plt.xlabel("Coefficient values")
            plt.ylabel("Frequency")
            fig.suptitle("Reliability coefficients histogram")
            fig.savefig(
                epi_fname[i].replace(
                    epi_fname[perm_volumes[i][0]].split("_")[-1],
                    epi_fname[perm_volumes[i][0]].split("_")[-1].split(".")[0]
                    + str(perm_volumes[i][0])
                    + str(perm_volumes[i][1])
                    + "_histogram.png",
                )
            )
            plt.close(fig)
            np.savetxt(
                    epi_fname[i].replace(
                        base_fname,
                        base_fname.split(".")[0]
                        + str(perm_volumes[i][0])
                        + str(perm_volumes[i][1])
                        + "_reliability.csv",
                    ),
                    reliability_dict[i],
                    delimiter=",",
                )

        if plot:
            plot_results = unmask(reliability_dict[i], mask)
            shape_epi = plot_results.shape
            sbref_epi = load_img(sbref)
            print("Loaded sbref for background: {sbref}")
            plot_results_affined = Nifti1Image(
                plot_results.get_fdata(),
                affine=sbref_epi.affine,
                header=sbref_epi.header,
            )
            print("Created new nilearn object to visualize results")
            title = (
                "Reliability map for "
                + epi_fname[perm_volumes[i][0]].split("_")[0].split("/")[-1]
                + " "
                + epi_fname[perm_volumes[i][0]].split("_")[-1].split(".")[0]
            )
            brain_reliability = plot_stat_map(
                plot_results_affined,
                sbref_epi,
                colorbar=True,
                draw_cross=False,
                title=title,
                cut_coords=(
                    (shape_epi[0] // 2),
                    (shape_epi[1] // 2),
                    (shape_epi[2] // 2),
                ),
                cmap="inferno",
                vmin=0,
                vmax=0.5,
            )
            brain_reliability.savefig(
                epi_fname[perm_volumes[i][0]].replace(
                    epi_fname[perm_volumes[i][0]].split("_")[-1],
                    epi_fname[perm_volumes[i][0]].split("_")[-1].split(".")[0]
                    + str(perm_volumes[i][0])
                    + str(perm_volumes[i][1])
                    + "_reliability.png",
                )
            )


#################
###### Main      ##############################
###############################################


source_dir = "/scratch/mflores/Resting_State/analysis_timeSeries"
methods = ["vanilla", "nordic", "tmmpca", "mppca", "nordic", "hydra"]
subjects = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-005"]
tasks = ["task-HABLA1200", "task-HABLA1700"]

for subject in subjects:
    for task in tasks:
        for method in methods:
            base_name = source_dir + "/" + subject + "_ses-1_" + task

            mask = base_name + "_echo-1_part-mag_gm_mask-union.nii.gz"

            sbref = (
                "/scratch/mflores/Resting_State/analysis/"
                + subject
                + "_ses-1_"
                + task
                + "_echo-1_part-mag_masked_sbref.nii.gz"
            )

            # EPI, Split
            try:
                epi = [base_name + "_OC_part-mag_bold_" + method + ".nii.gz"]
                print("LOG: Attempting EPI, split reliability analysis")
                reliability_analysis(epi, mask, sbref)
            except Exception:
                print(
                    f"ERROR: {subject}, task:{task}, and method:{method} one time series"
                )

            # EPI, Series
            try:
                epi_series = [
                    base_name + "_OC_part-mag_bold_" + method + "1.nii.gz",
                    base_name + "_OC_part-mag_bold_" + method + "2.nii.gz",
                ]

                print("LOG: Attempting EPI, series reliability analysis")
                reliability_analysis(epi_series, mask, sbref)
            except Exception:
                print(f"ERROR: {subject}, task:{task}, and method:{method} episeries")

            # Residuals, from split
            try:
                residual = [
                    base_name + "_residuals_part-mag_bold_" + method + ".nii.gz"
                ]
                print("LOG: Attempting residuals, split reliability analysis")
                reliability_analysis(residual, mask, sbref)
            except Exception:
                print(f"ERROR: {subject}, task:{task}, and method:{method} residual")

            # Residuals, from series
            try:
                residual_series = [
                    base_name + "_residuals_part-mag_bold_" + method + "1.nii.gz",
                    base_name + "_residuals_part-mag_bold_" + method + "2.nii.gz",
                ]
                print("LOG: Attempting residuals, series reliability analysis")
                reliability_analysis(residual_series, mask, sbref)
            except Exception:
                print("Error in residuals time series")
