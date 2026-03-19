#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:23:25 2023
@author: mflores

Compute voxel‑wise reliability for fMRI data.
Supports:
  - Single run: split into two halves.
  - Multiple runs: compare all unique pairs of runs.
Outputs:
  - Reliability NIfTI maps (optional).
  - Histograms of reliability values (optional).
  - CSV files of reliability values (optional).
"""

from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from nibabel import Nifti1Image
from nilearn.image import load_img
from nilearn.masking import apply_mask, unmask
from nilearn.plotting import plot_stat_map

# -----------------------------------------------------------
# Core function: voxel‑wise Pearson correlation (fast, vectorised)
# -----------------------------------------------------------
def voxelwise_correlation(X, Y):
    """
    Compute Pearson correlation per row (voxel) between two 2D arrays.
    Parameters
    ----------
    X, Y : 2D arrays of shape (n_voxels, n_timepoints)
    Returns
    -------
    r : 1D array of shape (n_voxels) with correlation coefficients.
    """
    # Z‑score along time axis
    Xz = (X - np.mean(X, axis=1, keepdims=True)) / np.std(X, axis=1, keepdims=True, ddof=1)
    Yz = (Y - np.mean(Y, axis=1, keepdims=True)) / np.std(Y, axis=1, keepdims=True, ddof=1)
    # Dot product and normalise by (n‑1)
    r = np.sum(Xz * Yz, axis=1) / (X.shape[1] - 1)
    return r

# -----------------------------------------------------------
# Main reliability analysis function
# -----------------------------------------------------------
def reliability_analysis(epi_fname, mask, sbref,
                         plot=False, savecorr=False, hist=False, make_nifti=True):
    """
    epi_fname : list of paths to NIfTI files (1 or more runs).
    mask      : path to binary mask NIfTI.
    sbref     : path to sbref image for background (used only if plot=True).
    plot      : if True, generate a stat map image of reliability.
    savecorr  : (kept for compatibility, but now saves reliability CSV instead).
    hist      : if True, plot and save histogram of reliability values.
    make_nifti: if True, save reliability map as NIfTI.
    """
    print("\n" + "="*50)
    print("Starting reliability analysis")
    print(f"EPI files: {epi_fname}")
    print(f"Mask: {mask}")
    print(f"Options: plot={plot}, savecorr={savecorr}, hist={hist}, make_nifti={make_nifti}")
    print("="*50)

    # -------------------------------------------------------
    # 1. Load and mask each run
    # -------------------------------------------------------
    array_dict = {}          # key = run index, value = masked data (voxels, time)
    for i, fname in enumerate(epi_fname):
        print(f"\nLoading {fname} ...")
        # apply_mask returns (time, voxels); we transpose to (voxels, time)
        data = np.transpose(apply_mask(fname, mask)).astype(np.float32)
        array_dict[i] = data
        print(f"   Shape (voxels, time): {data.shape}")
        n_voxels, n_time = data.shape
        print(f"   Voxels in mask: {n_voxels}")

    # Create submask to exclude voxels that are all zero (e.g., outside brain but inside mask)
    # Use first run as reference
    first_data = array_dict[0]
    array_submask = ~np.all(first_data == 0, axis=1)   # True for voxels with any non‑zero
    n_active = np.sum(array_submask)
    print(f"\nVoxels with non‑zero time series: {n_active}")

    # -------------------------------------------------------
    # 2. Determine run pairs and compute reliability per voxel
    # -------------------------------------------------------
    if len(array_dict) == 1:
        # Single run: split into two halves
        data = array_dict[0][array_submask, :]      # shape (n_active, time)
        mid = n_time // 2
        half1 = data[:, :mid]
        half2 = data[:, mid:]
        print(f"\nSplit single run at time point {mid} -> halves shape: {half1.shape}, {half2.shape}")
        runs = [0, 1]   # artificial run indices
        # Compute reliability between halves
        r_vals = voxelwise_correlation(half1, half2)
        reliability_dict = {0: r_vals}   # only one comparison
        pair_list = [(0, 1)]
        # For output filenames, we will use the original epi_fname[0] as base
        base_fname_for_output = epi_fname[0]

    else:
        # Multiple runs: compute reliability for every unique pair
        run_indices = list(array_dict.keys())
        pair_list = list(combinations(run_indices, 2))
        reliability_dict = {}
        print(f"\nFound {len(pair_list)} unique run pairs: {pair_list}")

        for idx, (i, j) in enumerate(pair_list):
            # Extract active voxels only
            data_i = array_dict[i][array_submask, :]
            data_j = array_dict[j][array_submask, :]
            print(f"  Computing pair {i}-{j}: shapes {data_i.shape}, {data_j.shape}")
            r_vals = voxelwise_correlation(data_i, data_j)
            reliability_dict[idx] = r_vals
        # For output filenames, we will use the first file of each pair as base
        # We'll handle this inside the loop

    # -------------------------------------------------------
    # 3. Generate outputs for each reliability map
    # -------------------------------------------------------
    for pair_idx, (i, j) in enumerate(pair_list if len(array_dict)>1 else [(0,1)]):
        r_vals = reliability_dict[pair_idx]   # 1D array of length n_active

        # ---- Save reliability values as CSV (if savecorr requested) ----
        if savecorr:
            csv_fname = epi_fname[i].replace(
                epi_fname[i].split('_')[-1],
                f"reliability_run{i}_run{j}.csv"
            )
            np.savetxt(csv_fname, r_vals, delimiter=',')
            print(f"Saved reliability CSV: {csv_fname}")

        # ---- Create NIfTI map (if make_nifti requested) ----
        if make_nifti:
            # Start with zeros for all voxels in mask
            vol_data = np.zeros(n_voxels, dtype=np.float32)
            vol_data[array_submask] = r_vals
            reliability_img = unmask(vol_data, mask)
            out_nii = epi_fname[i].replace(
                epi_fname[i].split('_')[-1],
                f"reliability_run{i}_run{j}.nii.gz"
            )
            reliability_img.to_filename(out_nii)
            print(f"Saved reliability NIfTI: {out_nii}")

        # ---- Histogram (if hist requested) ----
        if hist:
            fig, ax = plt.subplots(1, 1, figsize=(8,5))
            ax.hist(r_vals, bins=100, density=True, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Reliability (Pearson r)')
            ax.set_ylabel('Density')
            ax.set_title(f'Reliability histogram (runs {i}-{j})')
            hist_fname = epi_fname[i].replace(
                epi_fname[i].split('_')[-1],
                f"reliability_hist_run{i}_run{j}.png"
            )
            fig.savefig(hist_fname, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved histogram: {hist_fname}")

        # ---- Statistical map plot (if plot requested) ----
        if plot:
            # We need a volume for plotting; if make_nifti already created it, reuse.
            if not make_nifti:
                # Create temporary volume
                vol_data = np.zeros(n_voxels, dtype=np.float32)
                vol_data[array_submask] = r_vals
                reliability_img = unmask(vol_data, mask)

            # Load sbref for background
            sbref_img = load_img(sbref)
            # Ensure reliability image has same affine/header as sbref (for overlay)
            # (unmask already gives image with original affine, but sbref may differ)
            # Use sbref's affine for plotting background, but reliability data is already in its own space.
            # We'll plot reliability over sbref; if spaces differ, consider resampling (optional).
            # Here we simply use the reliability image as is, and provide sbref as background.
            display = plot_stat_map(
                reliability_img,
                bg_img=sbref_img,
                colorbar=True,
                draw_cross=False,
                title=f'Reliability (runs {i}-{j})',
                cut_coords=(0, 0, 0),   # adjust as needed
                cmap='inferno',
                vmin=0,
                vmax=0.5
            )
            plot_fname = epi_fname[i].replace(
                epi_fname[i].split('_')[-1],
                f"reliability_map_run{i}_run{j}.png"
            )
            display.savefig(plot_fname)
            display.close()
            print(f"Saved reliability map plot: {plot_fname}")

    print("\nReliability analysis completed.\n")


# -----------------------------------------------------------
# Main script (as provided by user)
# -----------------------------------------------------------
if __name__ == "__main__":
    source_dir = "/scratch/mflores/Resting_State/analysis_timeSeries"
    methods = ["vanilla", "nordic", "tmppca", "mppca", "nordic", "hydra"]
    subjects = ["sub-001", "sub-002", "sub-003", "sub-004", "sub-005"]
    tasks = ["task-HABLA1200", "task-HABLA1700"]
    runs = [""]   # no run label in filenames (adjust if needed)

    for subject in subjects:
        for task in tasks:
            for method in methods:
                for run in runs:
                    base_name = f"{source_dir}/{subject}_ses-1_{task}"
                    mask = f"{base_name}_echo-1_part-mag_gm_mask-union.nii.gz"
                    sbref = (f"/scratch/mflores/Resting_State/analysis/"
                             f"{subject}_ses-1_{task}_echo-1_part-mag_masked_sbref.nii.gz")

                    # ---- Single run, split analysis ----
                    try:
                        epi_single = [f"{base_name}_OC_part-mag_bold_{method}.nii.gz"]
                        print(f"\nLOG: Attempting EPI split analysis for {subject} {task} {method}")
                        reliability_analysis(epi_single, mask, sbref,
                                             plot=False, savecorr=False, hist=False, make_nifti=True)
                    except Exception as e:
                        print(f"ERROR (split) {subject}, {task}, {method}: {e}")

                    # ---- Two‑run series analysis ----
                    try:
                        epi_series = [
                            f"{base_name}_OC_part-mag_bold_{method}1.nii.gz",
                            f"{base_name}_OC_part-mag_bold_{method}2.nii.gz"
                        ]
                        print(f"\nLOG: Attempting EPI series analysis for {subject} {task} {method}")
                        reliability_analysis(epi_series, mask, sbref,
                                             plot=False, savecorr=False, hist=False, make_nifti=True)
                    except Exception as e:
                        print(f"ERROR (series) {subject}, {task}, {method}: {e}")
