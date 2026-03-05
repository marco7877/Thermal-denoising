#!/usr/bin/env python

import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs, get_data, new_img_like, resample_to_img
from sklearn.model_selection import ShuffleSplit

# -------------------------------
# Utilities
# -------------------------------

def percent_signal_change(data, mask):
    """Apply percent signal change scaling to masked voxels."""
    mask_flat = mask.ravel()
    Y = data.reshape(-1, data.shape[-1])
    mean = Y[mask_flat].mean(axis=1, keepdims=True)
    mean[mean == 0] = 1
    Y[mask_flat] = (Y[mask_flat] - mean) / mean * 100
    return Y.reshape(data.shape)

def zscore(X):
    """Column-wise z-score, avoid division by zero."""
    return (X - X.mean(axis=0)) / np.where(X.std(axis=0) == 0, 1, X.std(axis=0))

def block_diag_design(designs, n_task):
    """Create block-diagonal design matrix with run-specific nuisance regressors."""
    blocks = []
    for i, d in enumerate(designs):
        task = d.iloc[:, :n_task]
        nuisance = d.iloc[:, n_task:]
        nuisance = nuisance.add_suffix(f"_run{i}")
        blocks.append(pd.concat([task, nuisance], axis=1))
    return pd.concat(blocks).fillna(0).values

def compute_r2(Y_true, Y_pred):
    """Voxelwise R²."""
    ss_res = np.sum((Y_true - Y_pred) ** 2, axis=1)
    ss_tot = np.sum((Y_true - Y_true.mean(axis=1, keepdims=True)) ** 2, axis=1)
    return 1 - ss_res / (ss_tot + 1e-8)

def fast_glm(X, Y):
    """Solve OLS for all voxels at once."""
    XtX = X.T @ X
    XtY = X.T @ Y.T
    beta = np.linalg.solve(XtX, XtY)
    return beta.T

# -------------------------------
# Cross-validation GLM
# -------------------------------

def run_cv(nii_files, design_files, mask_img, n_task, tr, n_splits=50):

    mask = mask_img.get_fdata().astype(bool)

    # Load images
    imgs = [nib.load(f) for f in nii_files]
    designs = [pd.read_csv(f) for f in design_files]

    # Resample all images to reference affine if needed
    reference_img = imgs[0]
    imgs_resampled = [
        resample_to_img(img, reference_img, interpolation='continuous')
        if not np.allclose(img.affine, reference_img.affine, atol=1e-3) else img
        for img in imgs
    ]

    splitter = ShuffleSplit(n_splits=n_splits, test_size=0.3, random_state=0)

    r2_maps = []

    for train_idx, test_idx in splitter.split(imgs_resampled):

        # ---------- TRAIN ----------
        train_img = concat_imgs([imgs_resampled[i] for i in train_idx])
        train_data = percent_signal_change(get_data(train_img), mask)

        train_design = block_diag_design([designs[i] for i in train_idx], n_task)
        X_train = zscore(train_design)

        Y_train = train_data.reshape(-1, train_data.shape[-1])[mask.ravel()]
        Y_train = Y_train.T  # time × voxels

        beta = fast_glm(X_train, Y_train)
        task_betas = beta[:, :n_task]

        # ---------- TEST ----------
        test_img = concat_imgs([imgs_resampled[i] for i in test_idx])
        test_data = percent_signal_change(get_data(test_img), mask)

        test_design = block_diag_design([designs[i] for i in test_idx], n_task)
        X_test = zscore(test_design)[:, :n_task]

        Y_true = test_data.reshape(-1, test_data.shape[-1])[mask.ravel()]
        Y_pred = task_betas @ X_test.T

        r2 = compute_r2(Y_true, Y_pred)

        r2_map = np.zeros(mask.size)
        r2_map[mask.ravel()] = r2
        r2_maps.append(r2_map.reshape(mask.shape))

    return np.mean(r2_maps, axis=0)

# -------------------------------
# CLI
# -------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Fast cross-validated GLM for fMRI")

    parser.add_argument("--nii_files", nargs="+", required=True, help="List of NIfTI runs")
    parser.add_argument("--design_files", nargs="+", required=True, help="CSV design matrices per run")
    parser.add_argument("--mask", required=True, help="Mask NIfTI file")
    parser.add_argument("--n_task_regressors", type=int, required=True, help="Number of task regressors")
    parser.add_argument("--tr", type=float, default=2.0, help="Repetition time")
    parser.add_argument("--output", required=True, help="Output R² NIfTI filename")

    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    r2 = run_cv(args.nii_files, args.design_files, mask_img, args.n_task_regressors, args.tr)

    nib.save(nib.Nifti1Image(r2.astype(np.float32), mask_img.affine), args.output)
    print(f"Saved R² map to {args.output}")
