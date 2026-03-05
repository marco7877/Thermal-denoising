#!/usr/bin/env python
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs, get_data, resample_to_img
from sklearn.model_selection import ShuffleSplit
import argparse
import gc

# -------------------------------
# Utilities
# -------------------------------
def percent_signal_change(data, mask):
    """Apply percent signal change (PSC) scaling per voxel."""
    mask_flat = mask.ravel()
    Y = data.reshape(-1, data.shape[-1])
    voxel_means = Y[mask_flat].mean(axis=1, keepdims=True)
    voxel_means[voxel_means == 0] = 1
    Y[mask_flat] = (Y[mask_flat] - voxel_means) / voxel_means * 100
    return Y.reshape(data.shape)

def zscore(X):
    """Z-score design matrix (columns)."""
    return (X - X.mean(axis=0)) / np.where(X.std(axis=0) == 0, 1, X.std(axis=0))

def block_diag_design(designs, n_task):
    """
    Build a block-diagonal design matrix:
    - Task regressors shared across runs
    - Nuisance regressors run-specific
    """
    blocks = []
    for i, df in enumerate(designs):
        task = df.iloc[:, :n_task]
        nuisance = df.iloc[:, n_task:].add_suffix(f"_run{i}")
        blocks.append(pd.concat([task, nuisance], axis=1))
    return pd.concat(blocks).fillna(0).values.astype(np.float32)

def compute_r2(Y_true, Y_pred):
    """Compute R² per voxel."""
    ss_res = np.sum((Y_true - Y_pred)**2, axis=0)
    ss_tot = np.sum((Y_true - Y_true.mean(axis=0))**2, axis=0)
    return 1 - ss_res / (ss_tot + 1e-8)

# -------------------------------
# Cross-validation
# -------------------------------
def run_cv(
    nii_files,
    events_files,
    mask_img,
    n_task_regressors,
    tr,
    n_splits=50,
    test_size=0.3,
    random_state=42,
    permute=False
):
    mask = mask_img.get_fdata().astype(bool)
    mask_flat = mask.ravel()

    # Load and resample images
    imgs = [nib.load(f) for f in nii_files]
    reference_img = imgs[0]
    imgs_resampled = [
        resample_to_img(img, reference_img, interpolation='continuous')
        if not np.allclose(img.affine, reference_img.affine, atol=1e-3) else img
        for img in imgs
    ]

    # Load design matrices
    designs = [pd.read_csv(f).loc[:, ~pd.read_csv(f).columns.str.contains('^Unnamed')] for f in events_files]

    splitter = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)
    r2_maps = []

    for split_id, (train_idx, test_idx) in enumerate(splitter.split(imgs_resampled)):
        print(f"\n--- Split {split_id+1}/{n_splits} ---")

        # ----- TRAIN -----
        train_imgs = [imgs_resampled[i] for i in train_idx]
        train_img_concat = concat_imgs(train_imgs)
        train_data = percent_signal_change(get_data(train_img_concat), mask)
        X_train = block_diag_design([designs[i] for i in train_idx], n_task_regressors)
        Y_train = train_data[mask].T  # timepoints × voxels

        if permute:
            rng = np.random.default_rng(random_state + split_id)
            X_train[:,:n_task_regressors] = rng.permutation(X_train[:,:n_task_regressors])

        # Solve OLS
        beta = np.linalg.lstsq(X_train, Y_train, rcond=None)[0]  # regressors × voxels

        # ----- TEST -----
        test_imgs = [imgs_resampled[i] for i in test_idx]
        test_img_concat = concat_imgs(test_imgs)
        test_data = percent_signal_change(get_data(test_img_concat), mask)
        X_test = block_diag_design([designs[i] for i in test_idx], n_task_regressors)
        Y_true = test_data[mask]  # timepoints × voxels

        # Prediction (task regressors only)
        task_beta = beta[:n_task_regressors, :]
        Y_pred = X_test[:, :n_task_regressors] @ task_beta  # timepoints × voxels

        # R² per voxel
        r2 = compute_r2(Y_true, Y_pred)
        r2_map = np.zeros(mask.size)
        r2_map[mask_flat] = r2
        r2_maps.append(r2_map.reshape(mask.shape))

        # Cleanup
        del train_data, test_data, beta, Y_train, Y_true, Y_pred
        gc.collect()

    return np.mean(r2_maps, axis=0)

# -------------------------------
# CLI
# -------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_task_regressors", type=int, required=True)
    parser.add_argument("--tr", type=float, default=2.0)
    parser.add_argument("--n_splits", type=int, default=50)
    parser.add_argument("--test_size", type=float, default=0.3)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--permute", action="store_true")
    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    mean_r2 = run_cv(
        args.nii_files,
        args.events_files,
        mask_img,
        args.n_task_regressors,
        args.tr,
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state,
        permute=args.permute
    )

    out_file = f"{args.output_prefix}_mean_r2.nii.gz"
    nib.save(nib.Nifti1Image(mean_r2.astype(np.float32), mask_img.affine), out_file)
    print(f"Saved: {out_file}")
