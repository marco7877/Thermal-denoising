#!/usr/bin/env python3

import argparse
import numpy as np
import nibabel as nib
import pandas as pd

from nilearn.masking import apply_mask, unmask
from sklearn.model_selection import ShuffleSplit


# --------------------------------------------------
# Percent signal change (per run)
# --------------------------------------------------

def percent_signal_change(Y):

    mean = Y.mean(axis=1, keepdims=True)
    mean[mean == 0] = 1

    return (Y - mean) / mean * 100


# --------------------------------------------------
# Load fMRI runs
# --------------------------------------------------

def load_runs(nii_files, mask_img):

    runs = []

    for f in nii_files:

        img = nib.load(f)

        Y = apply_mask(img, mask_img).T   # voxels × time
        Y = percent_signal_change(Y)

        runs.append(Y)

    return runs


# --------------------------------------------------
# Load design matrices
# --------------------------------------------------

def load_designs(event_files):

    designs = []

    for f in event_files:

        df = pd.read_csv(f)
        df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

        designs.append(df.values.astype(np.float32))

    return designs


# --------------------------------------------------
# Build block design matrix
# --------------------------------------------------

def build_train_design(designs, train_idx, n_task):

    X_blocks = []

    n_runs = len(train_idx)

    nuis_cols = designs[train_idx[0]].shape[1] - n_task

    for r_i, r in enumerate(train_idx):

        X = designs[r]

        task = X[:, :n_task]
        nuis = X[:, n_task:]

        T = task.shape[0]

        nuis_block = np.zeros((T, nuis_cols * n_runs))

        start = r_i * nuis_cols
        end = start + nuis_cols

        nuis_block[:, start:end] = nuis

        X_block = np.hstack([task, nuis_block])

        X_blocks.append(X_block)

    return np.vstack(X_blocks)


# --------------------------------------------------
# Fit GLM
# --------------------------------------------------

def fit_glm(X, Y):

    XtX = X.T @ X
    XtY = X.T @ Y.T

    betas = np.linalg.solve(XtX, XtY)

    return betas


# --------------------------------------------------
# Regress nuisance from test data
# --------------------------------------------------

def regress_nuisance(Y, X_nuis):

    if X_nuis.shape[1] == 0:
        return Y

    XtX = X_nuis.T @ X_nuis
    XtY = X_nuis.T @ Y.T

    beta = np.linalg.solve(XtX, XtY)

    nuisance_pred = (X_nuis @ beta).T

    return Y - nuisance_pred


# --------------------------------------------------
# Vectorized R²
# --------------------------------------------------

def compute_r2(Y_true, Y_pred):

    ss_res = np.sum((Y_true - Y_pred) ** 2, axis=1)
    ss_tot = np.sum((Y_true - Y_true.mean(axis=1, keepdims=True)) ** 2, axis=1)

    r2 = 1 - ss_res / (ss_tot + 1e-8)

    return r2


# --------------------------------------------------
# Cross validation
# --------------------------------------------------

def run_cv(runs, designs, mask_img, splits, n_task):

    n_vox = runs[0].shape[0]

    r2_accum = np.zeros((len(splits), n_vox))

    for s_i, (train_idx, test_idx) in enumerate(splits):

        print(f"Split {s_i+1}/{len(splits)}")

        # -----------------------
        # training data
        # -----------------------

        X_train = build_train_design(designs, train_idx, n_task)

        Y_train = np.hstack([runs[i] for i in train_idx])

        betas = fit_glm(X_train, Y_train)

        task_betas = betas[:n_task].T

        # -----------------------
        # testing
        # -----------------------

        X_test_task = []
        X_test_nuis = []
        Y_test = []

        for r in test_idx:

            X = designs[r]

            task = X[:, :n_task]
            nuis = X[:, n_task:]

            X_test_task.append(task)
            X_test_nuis.append(nuis)

            Y_test.append(runs[r])

        X_test_task = np.vstack(X_test_task)
        X_test_nuis = np.vstack(X_test_nuis)

        Y_test = np.hstack(Y_test)

        # remove nuisance
        Y_clean = regress_nuisance(Y_test, X_test_nuis)

        # predict task
        Y_pred = (X_test_task @ task_betas).T

        r2 = compute_r2(Y_clean, Y_pred)

        r2_accum[s_i] = r2

    mean_r2 = r2_accum.mean(axis=0)

    return unmask(mean_r2, mask_img)


# --------------------------------------------------
# CLI
# --------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_task_regressors", type=int, required=True)

    parser.add_argument("--n_splits", type=int, default=20)
    parser.add_argument("--test_size", type=float, default=0.3)

    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    runs = load_runs(args.nii_files, mask_img)
    designs = load_designs(args.events_files)

    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=0
    )

    splits = list(splitter.split(runs))

    r2_img = run_cv(
        runs,
        designs,
        mask_img,
        splits,
        args.n_task_regressors
    )

    nib.save(r2_img, args.output)

    print("Saved:", args.output)


if __name__ == "__main__":
    main()