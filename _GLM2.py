#!/usr/bin/env python

import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs, get_data
from sklearn.model_selection import ShuffleSplit


# --------------------------------------------------
# Utilities
# --------------------------------------------------

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


# --------------------------------------------------
# Fast GLM solver
# --------------------------------------------------

def fast_glm(X, Y):

    XtX = X.T @ X
    XtY = X.T @ Y

    beta = np.linalg.solve(XtX, XtY)

    return beta.T


# --------------------------------------------------
# Cross validation
# --------------------------------------------------

def run_cv(nii_files, design_files, mask_img, n_task, tr, n_splits=50):

    mask = mask_img.get_fdata().astype(bool)

    imgs = [nib.load(f) for f in nii_files]
    designs = [pd.read_csv(f) for f in design_files]

    splitter = ShuffleSplit(n_splits=n_splits, test_size=0.3, random_state=0)

    r2_maps = []
    split_means = []

    for split_id, (train_idx, test_idx) in enumerate(splitter.split(imgs)):

        print(f"\n------ Split {split_id + 1}/{n_splits} ------")

        # ---------- TRAIN ----------
        train_img = concat_imgs([imgs[i] for i in train_idx])
        train_data = percent_signal_change(get_data(train_img), mask)

        train_design = block_diag_design(
            [designs[i] for i in train_idx],
            n_task
        )

        X_train = zscore(train_design)

        Y_train = train_data.reshape(-1, train_data.shape[-1])[mask.ravel()]
        Y_train = Y_train.T  # time × voxels

        beta = fast_glm(X_train, Y_train)

        task_betas = beta[:, :n_task]

        # ---------- TEST ----------
        test_img = concat_imgs([imgs[i] for i in test_idx])
        test_data = percent_signal_change(get_data(test_img), mask)

        test_design = block_diag_design(
            [designs[i] for i in test_idx],
            n_task
        )

        X_test_full = zscore(test_design)

        X_test = X_test_full[:, :n_task]
        X_nuisance = X_test_full[:, n_task:]

        Y_true = test_data.reshape(-1, test_data.shape[-1])[mask.ravel()]
        Y_true = Y_true.T  # time × voxels

        # --------------------------------------------------
        # Regress nuisance out of test data
        # --------------------------------------------------

        if X_nuisance.shape[1] > 0:
            beta_nuis = np.linalg.lstsq(X_nuisance, Y_true, rcond=None)[0]
            Y_true = Y_true - X_nuisance @ beta_nuis

        Y_true = Y_true.T  # voxels × time

        # --------------------------------------------------

        Y_pred = task_betas @ X_test.T

        # ---------- diagnostics ----------
        print("Y_true range:", np.min(Y_true), "to", np.max(Y_true))
        print("Y_pred range:", np.min(Y_pred), "to", np.max(Y_pred))

        r2 = compute_r2(Y_true, Y_pred)

        mean_r2 = np.mean(r2)
        split_means.append(mean_r2)

        print("Mean R2 this split:", mean_r2)
        print("R2 range:", np.min(r2), "to", np.max(r2))

        # ---------- map ----------
        r2_map = np.zeros(mask.size)
        r2_map[mask.ravel()] = r2
        r2_maps.append(r2_map.reshape(mask.shape))

    print("\n==============================")
    print("Cross-validation finished")
    print("Mean R2 across splits:", np.mean(split_means))
    print("Std R2 across splits:", np.std(split_means))
    print("==============================\n")

    return np.mean(r2_maps, axis=0)


# --------------------------------------------------
# CLI
# --------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--design_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_task_regressors", type=int, required=True)
    parser.add_argument("--tr", type=float, default=2.0)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    r2 = run_cv(
        args.nii_files,
        args.design_files,
        mask_img,
        args.n_task_regressors,
        args.tr
    )

    nib.save(
        nib.Nifti1Image(r2.astype(np.float32), mask_img.affine),
        args.output
    )
