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
def reliability_simulated(
        simulations_number, epi_fname, mask, sbref, plot=True, savecorr=False, hist=True, make_nifti=True, mu=0,sigma=1,seed=7877
        ):

    # Define variables
    array_dict = {} 
    corr_dict = {}
    reliability_dict = {}
    np.random.seed(seed)
    rng = np.random.default_rng()
    print("Worth double checking! To understand output ")
    print(
            f"Saving ... correlation matrixes: {savecorr}, r-values histogram: {hist}, plot: {plot}"
            )

    ##############################
    print("Getting timeseries shape")
    ##############################

    print(f"Loading epi file: {epi_fname} while applying mask: {mask}")
    array_dict[0] = np.transpose(apply_mask(epi_fname, mask))
    print(" Data loaded and masked!")
    shape = array_dict[0].shape
    print(f"Mask: {mask} contains {shape[0]} voxels")
    base_fname = epi_fname.split("_")[-1]

    # Simulation number random runs, compute correlation for each
    # Data is a  half normal between 0, and 1 with the shape of
    # epi_fname

    for i in range(simulations_number):
        corr_dict[i] = np.corrcoef(abs(rng.normal(mu,sigma,shape)))
        print(
                f"Functional connectivity for computed (pearson correlation) with shape {corr_dict[i].shape}"

                )
            # TODO should this be corr_dict[i]
        if savecorr:
            np.savetxt(
                    epi_fname.replace(
                        base_fname,
                        base_fname.split(".")[0]
                        + str(i)
                        + "simulated_connectivity.csv",
                        ),
                    corr_dict[i],
                    delimiter=",",
                        )

    perm_volumes = list(combinations(range(len(corr_dict)), 2)) 
    print(f" Calculating reliability for {len(perm_volumes)} combinations")

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
                        + str(perm_volumes[i][1])
                        + "simulated",
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
                        + "simulated_histogram.png",
                        )
                    )
            plt.close(fig)
            np.savetxt(
                    epi_fname[i].replace(
                        base_fname,
                        base_fname.split(".")[0]
                        + str(perm_volumes[i][0])
                        + str(perm_volumes[i][1])
                        + "simulated_reliability.csv",
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
                    "Reliability map for simulated data "
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

method = "vanilla"

subject = "sub-001"

task = "task-HABLA1200"

base_name = source_dir + "/" + subject + "_ses-1_" + task

mask = base_name + "_echo-1_part-mag_gm_mask-union.nii.gz"

sbref = (
        "/scratch/mflores/Resting_State/analysis/"
        + subject
        + "_ses-1_"
        + task
        + "_echo-1_part-mag_masked_sbref.nii.gz"
        )

            # EPI, Simulate
try:
    epi = [base_name + "_OC_part-mag_bold_" + method + ".nii.gz"]
    print("LOG: Attempting SIMULATION, reliability analysis")
    reliability_simulated(100, epi, mask, sbref)
except Exception:
        print(
                f"ERROR: {subject}, task:{task}, and method:{method} one time series"
                )
