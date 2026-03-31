#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from nibabel import load, Nifti1Image
from nilearn.image import concat_imgs, get_data, new_img_like, resample_to_img
from nilearn.masking import apply_mask
from nilearn.glm.first_level import FirstLevelModel
import warnings
warnings.filterwarnings('ignore')


def percent_change_scaling_per_run(img, mask_img, debug=False):
    """
    Apply percent change scaling PER RUN: (x - mean)/mean * 100
    This is done independently for each run before concatenation.
    """
    data = get_data(img).astype(np.float64)
    mask_data = mask_img.get_fdata().astype(bool)

    original_shape = data.shape
    n_timepoints = original_shape[-1]
    n_voxels = np.prod(original_shape[:-1])
    data_reshaped = data.reshape(n_voxels, n_timepoints)
    mask_flat = mask_data.ravel()

    masked_data = data_reshaped[mask_flat, :].copy()
    if masked_data.size > 0:
        voxel_means = np.mean(masked_data, axis=1, keepdims=True)
        voxel_means_safe = np.where(np.abs(voxel_means) < 1e-6, 1, voxel_means)
        masked_pc = ((masked_data - voxel_means_safe) / voxel_means_safe) * 100
        data_pc_reshaped = np.zeros_like(data_reshaped)
        data_pc_reshaped[mask_flat, :] = masked_pc
    else:
        data_pc_reshaped = data_reshaped

    data_pc = data_pc_reshaped.reshape(original_shape).astype(np.float32)
    return new_img_like(img, data_pc)


def load_and_preprocess_runs(nii_files, mask_img, debug=False):
    """Load all runs, optionally resample to first run, and apply percent change scaling."""
    all_runs_imgs = []
    run_timepoints = []
    reference_img = None

    for f in nii_files:
        img = load(f)
        if reference_img is None:
            reference_img = img
        else:
            if not np.allclose(img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
                img = resample_to_img(img, reference_img, interpolation='continuous')
        run_timepoints.append(img.shape[-1])
        img_pc = percent_change_scaling_per_run(img, mask_img, debug=debug)
        all_runs_imgs.append(img_pc)

    return all_runs_imgs, reference_img, run_timepoints


def load_design_matrices(events_files, run_timepoints, n_task_regressors):
    """Load precomputed design matrices and verify they match run timepoints."""
    designs = []
    for i, f in enumerate(events_files):
        df = pd.read_csv(f)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]  # drop empty columns
        if len(df) != run_timepoints[i]:
            raise ValueError(f"Timepoint mismatch for run {i}: fMRI={run_timepoints[i]}, Design={len(df)}")
        designs.append(df)
    return designs


def build_block_design(designs, n_task_regressors):
    """
    Build a block‑diagonal design matrix from a list of run‑specific DataFrames.
    Task columns are shared across runs; nuisance columns are renamed per run.
    """
    blocks = []
    for run_idx, df in enumerate(designs):
        task_part = df.iloc[:, :n_task_regressors]
        nuisance_part = df.iloc[:, n_task_regressors:]
        # Rename nuisance columns to be run‑specific
        rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
        nuisance_part = nuisance_part.rename(columns=rename_dict)
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        blocks.append(run_design)
    return pd.concat(blocks, axis=0, ignore_index=True).fillna(0)


def main():
    parser = argparse.ArgumentParser(description="Simple first‑level GLM")
    parser.add_argument("--nii_files", nargs="+", required=True, help="List of NIfTI files (one per run)")
    parser.add_argument("--events_files", nargs="+", required=True, help="List of CSV design matrices (one per run)")
    parser.add_argument("--mask", required=True, help="Mask image")
    parser.add_argument("--n_task_regressors", type=int, required=True,
                        help="Number of task regressors at the beginning of each design matrix")
    parser.add_argument("--contrast", type=str,
                        default="face1-(bodylimb1+ES_SC1+ES_RW1+ES_FF1+ES_CS1)/5",
                        help="Contrast expression (use column names from the design matrix)")
    parser.add_argument("--tr", type=float, default=2.0, help="Repetition time (seconds)")
    parser.add_argument("--hrf_model", type=str, default="spm", help="HRF model (e.g., 'spm', 'glover')")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--debug", action="store_true", help="Print debug information")
    args = parser.parse_args()

    # Validate inputs
    if len(args.nii_files) != len(args.events_files):
        raise ValueError("Number of NIfTI files and event files must match")

    # Load mask
    mask_img = load(args.mask)
    mask_data = mask_img.get_fdata().astype(bool)

    # Load and preprocess runs (including per‑run percent change scaling)
    all_runs_imgs, _, run_timepoints = load_and_preprocess_runs(
        args.nii_files, mask_img, debug=args.debug
    )

    # Load design matrices
    designs = load_design_matrices(args.events_files, run_timepoints, args.n_task_regressors)

    # Build block‑diagonal design matrix
    design_matrix = build_block_design(designs, args.n_task_regressors)

    # Concatenate all runs into one 4D image
    concat_img = concat_imgs(all_runs_imgs)

    # Fit GLM
    print("Fitting GLM...")
    glm = FirstLevelModel(
        t_r=args.tr,
        mask_img=None,
        standardize=True,
        signal_scaling=False,
        hrf_model=args.hrf_model,
        minimize_memory=True,
    )
    glm = glm.fit(concat_img, design_matrices=design_matrix)

    # Compute contrast
    print(f"Computing contrast: {args.contrast}")
    t_img = glm.compute_contrast(args.contrast, stat_type="t")
    beta_img = glm.compute_contrast(args.contrast, stat_type="effect_size")   # effect size (beta)

    # Save outputs
    t_out = f"{args.output_prefix}_t.nii.gz"
    beta_out = f"{args.output_prefix}_beta.nii.gz"
    t_img.to_filename(t_out)
    beta_img.to_filename(beta_out)
    print(f"Saved t‑statistic map: {t_out}")
    print(f"Saved beta map: {beta_out}")

    # Optional: also save the design matrix used (for verification)
    # design_matrix.to_csv(f"{args.output_prefix}_design_matrix.csv", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
