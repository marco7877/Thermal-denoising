#!/usr/bin/env python

import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs, get_data, resample_to_img
from sklearn.model_selection import ShuffleSplit

# -------------------------------
# Utilities
# -------------------------------

def percent_signal_change(data, mask):
    mask_flat = mask.ravel()
    Y = data.reshape(-1, data.shape[-1])
    mean = Y[mask_flat].mean(axis=1, keepdims=True)
    mean[mean == 0] = 1
    Y[mask_flat] = (Y[mask_flat] - mean) / mean * 100
    return Y.reshape(data.shape)

def zscore(X):
    return (X - X.mean(axis=0)) / np.where(X.std(axis=0) == 0, 1, X.std(axis=0))

def block_diag_design(designs, n_task):
    blocks = []
    for i, d in enumerate(designs):
        task = d.iloc[:, :n_task]
        nuisance = d.iloc[:, n_task:]
        nuisance = nuisance.add_suffix(f"_run{i}")
        blocks.append(pd.concat([task, nuisance], axis=1))
    return pd.concat(blocks).fillna(0).values

def compute_r2(Y_true, Y_pred):
    ss_res = np.sum((Y_true - Y_pred) ** 2, axis=1)
    ss_tot = np.sum((Y_true - Y_true.mean(axis=1, keepdims=True)) ** 2, axis=1)
    return 1 - ss_res / (ss_tot + 1e-8)

def fast_glm(X, Y):
    XtX = X.T @ X
    XtY = X.T @ Y.T
    beta = np.linalg.solve(XtX, XtY)
    return beta.T

# -------------------------------
# Cross-validation GLM
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

    imgs = [nib.load(f) for f in nii_files]
    designs = [pd.read_csv(f) for f in events_files]

    # Resample all images to reference
    reference_img = imgs[0]
    imgs_resampled = [
        resample_to_img(img, reference_img, interpolation='continuous')
        if not np.allclose(img.affine, reference_img.affine, atol=1e-3) else img
        for img in imgs
    ]

    splitter = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)

    r2_maps = []

    rng = np.random.default_rng(random_state)

    for train_idx, test_idx in splitter.split(imgs_resampled):

        # ---------- TRAIN ----------
        train_img = concat_imgs([imgs_resampled[i] for i in train_idx])
        train_data = percent_signal_change(get_data(train_img), mask)

        train_design_list = [designs[i] for i in train_idx]
        if permute:
            train_design_list = [d.sample(frac=1, random_state=rng).reset_index(drop=True) for d in train_design_list]

        train_design = block_diag_design(train_design_list, n_task_regressors)
        X_train = zscore(train_design)

        Y_train = train_data.reshape(-1, train_data.shape[-1])[mask.ravel()]
        Y_train = Y_train.T  # time × voxels

        beta = fast_glm(X_train, Y_train)
        task_betas = beta[:, :n_task_regressors]

        # ---------- TEST ----------
        test_img = concat_imgs([imgs_resampled[i] for i in test_idx])
        test_data = percent_signal_change(get_data(test_img), mask)

        test_design_list = [designs[i] for i in test_idx]
        if permute:
            test_design_list = [d.sample(frac=1, random_state=rng).reset_index(drop=True) for d in test_design_list]

        test_design = block_diag_design(test_design_list, n_task_regressors)
        X_test = zscore(test_design)[:, :n_task_regressors]

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

    parser = argparse.ArgumentParser(description="Cross-validated GLM with block-diagonal design matrices")

    parser.add_argument("--nii_files", nargs="+", required=True, help="List of NIfTI runs")
    parser.add_argument("--events_files", nargs="+", required=True, help="CSV design matrices per run")
    parser.add_argument("--mask", required=True, help="Mask NIfTI file")
    parser.add_argument("--n_task_regressors", type=int, required=True, help="Number of task regressors")
    parser.add_argument("--tr", type=float, default=2.0, help="Repetition time")
    parser.add_argument("--n_splits", type=int, default=50, help="Number of cross-validation splits")
    parser.add_argument("--test_size", type=float, default=0.3, help="Proportion of runs for testing")
    parser.add_argument("--n_permutations", type=int, default=0, help="Number of permutations for null distribution")
    parser.add_argument("--random_state", type=int, default=42, help="Random seed")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--only_permutations", action="store_true", help="Only run permutations")

    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    # Run real data CV
    if not args.only_permutations:
        mean_r2 = run_cv(
            args.nii_files,
            args.events_files,
            mask_img,
            args.n_task_regressors,
            args.tr,
            n_splits=args.n_splits,
            test_size=args.test_size,
            random_state=args.random_state,
            permute=False
        )

        out_mean = f"{args.output_prefix}_real_mean_r2.nii.gz"
        nib.save(nib.Nifti1Image(mean_r2.astype(np.float32), mask_img.affine), out_mean)
        print(f"Saved: {out_mean}")

    # Run permutations if requested
    if args.n_permutations > 0:
        max_distribution = np.zeros(args.n_permutations, dtype=np.float32)
        running_mean = None

        for p in range(args.n_permutations):
            mean_r2_perm = run_cv(
                args.nii_files,
                args.events_files,
                mask_img,
                args.n_task_regressors,
                args.tr,
                n_splits=args.n_splits,
                test_size=args.test_size,
                random_state=args.random_state + p + 1,
                permute=True
            )

            if running_mean is None:
                running_mean = np.zeros_like(mean_r2_perm, dtype=np.float32)

            running_mean += mean_r2_perm
            mask_data = mask_img.get_fdata().astype(bool)
            max_distribution[p] = np.max(mean_r2_perm[mask_data])

        perm_mean = running_mean / args.n_permutations
        out_perm = f"{args.output_prefix}_perm_mean_r2.nii.gz"
        nib.save(nib.Nifti1Image(perm_mean.astype(np.float32), mask_img.affine), out_perm)
        print(f"Saved: {out_perm}")

        out_dist = f"{args.output_prefix}_perm_max_distribution.npy"
        np.save(out_dist, max_distribution)
        print(f"Saved: {out_dist}")
