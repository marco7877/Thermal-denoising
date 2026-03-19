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
import os

# -----------------------------------------------------------
# Core function: voxel‑wise Pearson correlation (fast, vectorised)
# -----------------------------------------------------------
def voxelwise_correlation(X, Y):
    """
    Compute Pearson correlation per row (voxel) between two 2D arrays.

    Parameters
    ----------
    X, Y : 2D arrays of shape (n_voxels, n_timepoints)
           Must have the same number of time points.

    Returns
    -------
    r : 1D array of shape (n_voxels) with correlation coefficients.
    """
    # Z‑score along time axis (with ddof=1 for sample std)
    Xz = (X - np.mean(X, axis=1, keepdims=True)) / np.std(X, axis=1, keepdims=True, ddof=1)
    Yz = (Y - np.mean(Y, axis=1, keepdims=True)) / np.std(Y, axis=1, keepdims=True, ddof=1)
    # Dot product and normalise by (n‑1)
    r = np.sum(Xz * Yz, axis=1) / (X.shape[1] - 1)
    return r


def add_suffix_to_filename(fpath, suffix):
    """
    Insert a suffix before the file extension.
    Handles .nii.gz correctly.
    Example: add_suffix('sub-1_task-bold.nii.gz', '_reliability')
             -> 'sub-1_task-bold_reliability.nii.gz'
    """
    base, ext = os.path.splitext(fpath)
    if ext == '.gz':
        # .nii.gz case: base is e.g. 'sub-1_task-bold.nii'
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext   # now ext = '.nii.gz'
    return base + suffix + ext


def reliability_analysis(epi_fname, mask, sbref,
                         plot=False, savecorr=False, hist=False, make_nifti=True):
    """
    epi_fname : list of paths to NIfTI files (1 or more runs).
    mask      : path to binary mask NIfTI.
    sbref     : path to sbref image for background (used only if plot=True).
    plot      : if True, generate a stat map image of reliability.
    savecorr  : if True, save reliability values as CSV (one per voxel).
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
        print(f"   Voxels in mask: {data.shape[0]}")

    # Create submask to exclude voxels that are all zero (e.g., outside brain but inside mask)
    # Use first run as reference
    first_data = array_dict[0]
    array_submask = ~np.all(first_data == 0, axis=1)   # True for voxels with any non‑zero
    n_active = np.sum(array_submask)
    print(f"\nVoxels with non‑zero time series (first run): {n_active}")

    # -------------------------------------------------------
    # 2. Determine run pairs and compute reliability per voxel
    # -------------------------------------------------------
    if len(array_dict) == 1:
        # Single run: split into two halves
        data = array_dict[0][array_submask, :]      # shape (n_active, time)
        n_time_total = data.shape[1]
        mid = n_time_total // 2
        half1 = data[:, :mid]
        half2 = data[:, mid:]
        print(f"\nSplit single run at time point {mid} -> halves shape: {half1.shape}, {half2.shape}")
        # Halves may have slightly different lengths if odd number of volumes; truncate to min length
        min_len = min(half1.shape[1], half2.shape[1])
        if half1.shape[1] != half2.shape[1]:
            print(f"  Warning: halves have different lengths ({half1.shape[1]} vs {half2.shape[1]}). Truncating to {min_len}.")
            half1 = half1[:, :min_len]
            half2 = half2[:, :min_len]
        # Compute reliability between halves
        r_vals = voxelwise_correlation(half1, half2)
        reliability_dict = {0: r_vals}   # only one comparison
        pair_list = [(0, 1)]

    else:
        # Multiple runs: compute reliability for every unique pair
        run_indices = list(array_dict.keys())
        pair_list = list(combinations(run_indices, 2))
        reliability_dict = {}
        print(f"\nFound {len(pair_list)} unique run pairs: {pair_list}")

        for idx, (i, j) in enumerate(pair_list):
            # Extract active voxels only (based on first run's submask)
            data_i = array_dict[i][array_submask, :]
            data_j = array_dict[j][array_submask, :]
            print(f"  Computing pair {i}-{j}: shapes {data_i.shape}, {data_j.shape}")

            # Handle possible different lengths
            if data_i.shape[1] != data_j.shape[1]:
                min_len = min(data_i.shape[1], data_j.shape[1])
                print(f"    Warning: runs have different lengths. Truncating to {min_len} time points.")
                data_i = data_i[:, :min_len]
                data_j = data_j[:, :min_len]

            r_vals = voxelwise_correlation(data_i, data_j)
            reliability_dict[idx] = r_vals

    # -------------------------------------------------------
    # 3. Generate outputs for each reliability map
    # -------------------------------------------------------
    # Decide suffix format based on number of runs
    if len(array_dict) > 2:
        # For more than two runs, include pair indices to avoid overwriting
        suffix_template = '_reliability_run{}_run{}'
    else:
        # For one or two runs, just '_reliability'
        suffix_template = '_reliability'

    for pair_idx, (i, j) in enumerate(pair_list if len(array_dict)>1 else [(0,1)]):
        r_vals = reliability_dict[pair_idx]   # 1D array of length n_active

        # Base filename for outputs: use the first file of the pair
        base_fname = epi_fname[i]

        # Build the suffix for this pair
        if len(array_dict) > 2:
            suffix = suffix_template.format(i, j)
        else:
            suffix = suffix_template

        # ---- Save reliability values as CSV (if savecorr requested) ----
        if savecorr:
            csv_fname = add_suffix_to_filename(base_fname, suffix + '.csv')
            np.savetxt(csv_fname, r_vals, delimiter=',')
            print(f"Saved reliability CSV: {csv_fname}")

        # ---- Create NIfTI map (if make_nifti requested) ----
        if make_nifti:
            # Start with zeros for all voxels in mask
            vol_data = np.zeros(first_data.shape[0], dtype=np.float32)
            vol_data[array_submask] = r_vals
            reliability_img = unmask(vol_data, mask)
            out_nii = add_suffix_to_filename(base_fname, suffix + '.nii.gz')
            reliability_img.to_filename(out_nii)
            print(f"Saved reliability NIfTI: {out_nii}")

        # ---- Histogram (if hist requested) ----
        if hist:
            fig, ax = plt.subplots(1, 1, figsize=(8,5))
            ax.hist(r_vals, bins=100, density=True, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Reliability (Pearson r)')
            ax.set_ylabel('Density')
            ax.set_title(f'Reliability histogram (runs {i}-{j})')
            hist_fname = add_suffix_to_filename(base_fname, suffix + '_hist.png')
            fig.savefig(hist_fname, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved histogram: {hist_fname}")

        # ---- Statistical map plot (if plot requested) ----
        if plot:
            # We need a volume for plotting; if make_nifti already created it, reuse.
            if not make_nifti:
                # Create temporary volume
                vol_data = np.zeros(first_data.shape[0], dtype=np.float32)
                vol_data[array_submask] = r_vals
                reliability_img = unmask(vol_data, mask)

            # Load sbref for background
            sbref_img = load_img(sbref)
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
            plot_fname = add_suffix_to_filename(base_fname, suffix + '_map.png')
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