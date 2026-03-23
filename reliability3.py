#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute reliability of functional connectivity patterns for fMRI data,
using a batched approach to avoid memory blow‑up.
"""

from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from nibabel import Nifti1Image
from nilearn.image import load_img
from nilearn.masking import apply_mask, unmask
from nilearn.plotting import plot_stat_map
from scipy.stats import pearsonr
import os

# -----------------------------------------------------------
# Batched computation of reliability per voxel
# -----------------------------------------------------------
def reliability_from_connectivity_profiles_batched(X, Y, batch_size=1000):
    """
    X, Y : 2D arrays of shape (n_voxels, n_time) standardized (z‑scored)
           for run1 and run2 respectively.
    batch_size : number of voxels to process at once.

    Returns reliability array of length n_voxels (r²).
    """
    n_voxels, n_time = X.shape
    reliability = np.zeros(n_voxels, dtype=np.float32)

    denom = n_time - 1.0

    for start in range(0, n_voxels, batch_size):
        end = min(start + batch_size, n_voxels)
        batch_voxels = slice(start, end)

        # Correlation of batch voxels with all voxels in run1
        X_batch = X[batch_voxels, :]
        corr_batch1 = np.dot(X_batch, X.T) / denom   # shape (batch_size, n_voxels)

        Y_batch = Y[batch_voxels, :]
        corr_batch2 = np.dot(Y_batch, Y.T) / denom

        for idx_in_batch, global_idx in enumerate(range(start, end)):
            prof1 = corr_batch1[idx_in_batch, :]
            prof2 = corr_batch2[idx_in_batch, :]

            # Remove self‑correlation (which is 1.0)
            #prof1 = np.delete(prof1, global_idx)
            #prof2 = np.delete(prof2, global_idx)

            # Check for constant profile (should not happen due to mask, but safe)
            if np.std(prof1) == 0 or np.std(prof2) == 0:
                reliability[global_idx] = 0.0
                continue

            # Pearson correlation using scipy (could also compute directly)
            # Use only finite values
            mask = np.isfinite(prof1) & np.isfinite(prof2)
            if not np.all(mask):
                prof1 = prof1[mask]
                prof2 = prof2[mask]
            if len(prof1) < 2:
                reliability[global_idx] = 0.0
                continue

            r_val, _ = pearsonr(prof1, prof2)
            if np.isnan(r_val):
                reliability[global_idx] = 0.0
            else:
                r_val = np.clip(r_val, -1.0, 1.0)
                reliability[global_idx] = r_val ** 2

    return reliability


# -----------------------------------------------------------
# Helper functions
# -----------------------------------------------------------
def add_suffix_to_filename(fpath, suffix):
    """
    Insert a suffix before the file extension, handling .nii.gz.
    Example: add_suffix('sub-1_bold.nii.gz', '_reliability')
             -> 'sub-1_bold_reliability.nii.gz'
    """
    base, ext = os.path.splitext(fpath)
    if ext == '.gz':
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base + suffix + ext


def reliability_analysis(epi_fname, mask, sbref,
                         plot=False, savecorr=False, hist=False, make_nifti=True,
                         batch_size=1000):
    """
    Compute reliability of functional connectivity patterns using batched correlation.

    epi_fname : list of paths to NIfTI files (1 or more runs).
    mask      : path to binary mask NIfTI.
    sbref     : path to sbref image for background (used only if plot=True).
    plot      : if True, generate a stat map image of reliability.
    savecorr  : if True, save reliability values as CSV (one per voxel).
    hist      : if True, plot and save histogram of reliability values.
    make_nifti: if True, save reliability map as NIfTI.
    batch_size: number of voxels to process at once (adjust for memory).
    """
    print("\n" + "="*50)
    print("Starting reliability analysis (connectivity‑based, batched)")
    print(f"EPI files: {epi_fname}")
    print(f"Mask: {mask}")
    print("="*50)

    # -------------------------------------------------------
    # 1. Load and mask each run
    # -------------------------------------------------------
    array_dict = {}
    for i, fname in enumerate(epi_fname):
        print(f"\nLoading {fname} ...")
        data = np.transpose(apply_mask(fname, mask)).astype(np.float32)
        array_dict[i] = data
        print(f"   Shape (voxels, time): {data.shape}")

    # Active voxels (non‑zero in first run)
    first_data = array_dict[0]
    array_submask = ~np.all(first_data == 0, axis=1)
    n_active = np.sum(array_submask)
    print(f"\nVoxels with non‑zero time series: {n_active}")

    # -------------------------------------------------------
    # 2. Build run pairs and compute reliability per voxel
    # -------------------------------------------------------
    if len(array_dict) == 1:
        # Single run: split into two halves
        data = array_dict[0][array_submask, :]
        n_time = data.shape[1]
        mid = n_time // 2
        half1 = data[:, :mid]
        half2 = data[:, mid:]
        print(f"\nSplit into halves at time {mid}: shapes {half1.shape}, {half2.shape}")
        # Truncate to the shorter half
        min_len = min(half1.shape[1], half2.shape[1])
        if half1.shape[1] != half2.shape[1]:
            print(f"  Warning: halves have different lengths – truncating to {min_len}")
            half1 = half1[:, :min_len]
            half2 = half2[:, :min_len]

        # Standardize (z‑score) along time axis
        print("  Standardizing time series...")
        half1_z = (half1 - np.mean(half1, axis=1, keepdims=True)) / np.std(half1, axis=1, keepdims=True, ddof=1)
        half2_z = (half2 - np.mean(half2, axis=1, keepdims=True)) / np.std(half2, axis=1, keepdims=True, ddof=1)

        # Compute reliability in batches
        reliability = reliability_from_connectivity_profiles_batched(half1_z, half2_z, batch_size)
        reliability_dict = {0: reliability}
        pair_list = [(0, 1)]

    else:
        run_indices = list(array_dict.keys())
        pair_list = list(combinations(run_indices, 2))
        reliability_dict = {}
        print(f"\nFound {len(pair_list)} run pairs: {pair_list}")

        for idx, (i, j) in enumerate(pair_list):
            data_i = array_dict[i][array_submask, :]
            data_j = array_dict[j][array_submask, :]
            print(f"\n  Pair {i}-{j}: shapes {data_i.shape}, {data_j.shape}")

            # Ensure same number of time points
            if data_i.shape[1] != data_j.shape[1]:
                min_len = min(data_i.shape[1], data_j.shape[1])
                print(f"    Warning: runs have different lengths – truncating to {min_len}")
                data_i = data_i[:, :min_len]
                data_j = data_j[:, :min_len]

            # Standardize
            print("    Standardizing time series...")
            data_i_z = (data_i - np.mean(data_i, axis=1, keepdims=True)) / np.std(data_i, axis=1, keepdims=True, ddof=1)
            data_j_z = (data_j - np.mean(data_j, axis=1, keepdims=True)) / np.std(data_j, axis=1, keepdims=True, ddof=1)

            # Compute reliability in batches
            reliability = reliability_from_connectivity_profiles_batched(data_i_z, data_j_z, batch_size)
            reliability_dict[idx] = reliability

    # -------------------------------------------------------
    # 3. Generate outputs for each reliability map
    # -------------------------------------------------------
    if len(array_dict) > 2:
        suffix_template = '_reliability_conn_run{}_run{}'
    else:
        suffix_template = '_reliability_conn'

    for pair_idx, (i, j) in enumerate(pair_list if len(array_dict)>1 else [(0,1)]):
        reliability = reliability_dict[pair_idx]
        base_fname = epi_fname[i]

        if len(array_dict) > 2:
            suffix = suffix_template.format(i, j)
        else:
            suffix = suffix_template

        # Print summary of reliability values
        print(f"\nReliability for pair {i}-{j}: min={reliability.min():.6f}, max={reliability.max():.6f}, mean={reliability.mean():.6f}")

        # ---- Save CSV (if requested) ----
        if savecorr:
            base_without_ext, ext = os.path.splitext(base_fname)
            if ext == '.gz':
                base_without_ext, ext2 = os.path.splitext(base_without_ext)
            csv_fname = base_without_ext + suffix + '.csv'
            np.savetxt(csv_fname, reliability, delimiter=',')
            print(f"Saved reliability CSV: {csv_fname}")

        # ---- Save NIfTI map (if requested) ----
        if make_nifti:
            vol_data = np.zeros(first_data.shape[0], dtype=np.float32)
            vol_data[array_submask] = reliability
            reliability_img = unmask(vol_data, mask)
            out_nii = add_suffix_to_filename(base_fname, suffix)   # no extra extension
            reliability_img.to_filename(out_nii)
            print(f"Saved reliability NIfTI: {out_nii}")

        # ---- Histogram (if requested) ----
        if hist:
            base_without_ext, ext = os.path.splitext(base_fname)
            if ext == '.gz':
                base_without_ext, ext2 = os.path.splitext(base_without_ext)
            hist_fname = base_without_ext + suffix + '_hist.png'
            fig, ax = plt.subplots(1, 1, figsize=(8,5))
            ax.hist(reliability, bins=100, density=True, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Reliability (r²)')
            ax.set_ylabel('Density')
            ax.set_title(f'Reliability histogram (runs {i}-{j})')
            fig.savefig(hist_fname, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved histogram: {hist_fname}")

        # ---- Statistical map plot (if requested) ----
        if plot:
            if not make_nifti:
                vol_data = np.zeros(first_data.shape[0], dtype=np.float32)
                vol_data[array_submask] = reliability
                reliability_img = unmask(vol_data, mask)
            sbref_img = load_img(sbref)
            display = plot_stat_map(
                reliability_img,
                bg_img=sbref_img,
                colorbar=True,
                draw_cross=False,
                title=f'Connectivity reliability (runs {i}-{j})',
                cut_coords=(0, 0, 0),   # adjust as needed
                cmap='inferno',
                vmin=0,
                vmax=0.5
            )
            base_without_ext, ext = os.path.splitext(base_fname)
            if ext == '.gz':
                base_without_ext, ext2 = os.path.splitext(base_without_ext)
            plot_fname = base_without_ext + suffix + '_map.png'
            display.savefig(plot_fname)
            display.close()
            print(f"Saved reliability map plot: {plot_fname}")

    print("\nReliability analysis completed.\n")



# -----------------------------------------------------------
# Main script (as provided by user)
# -----------------------------------------------------------
if __name__ == "__main__":
    source_dir = "/scratch/mflores/Rest_HighRes/analysis_timeSeries"
#    methods = ["vanilla", "nordic", "tmppca", "mppca", "nordic", "hydra"]
    methods = ["mppca"]
    subjects = ["sub-001"]
    tasks = ["task-REST"]
    runs = ["_run-1"]   # no run label in filenames (adjust if needed)

    for subject in subjects:
        for task in tasks:
            for method in methods:
                for run in runs:
                    base_name = f"{source_dir}/{subject}_ses-1_{task}{run}"
                    mask = f"{base_name}_echo-1_part-mag_gm_mask-union.nii.gz"
                    sbref = (f"/scratch/mflores/Resting_State/analysis/"
                             f"{subject}_ses-1_{task}{run}_echo-1_part-mag_masked_sbref.nii.gz")

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
