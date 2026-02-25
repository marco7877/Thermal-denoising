#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from sklearn.model_selection import ShuffleSplit
from tqdm import tqdm
import gc


# ============================================================
# Utilities
# ============================================================

def compute_r2(y_true, y_pred):
    """
    Proper cross-validated R²
    Memory safe and float32
    """
    y_true = y_true.astype(np.float32)
    y_pred = y_pred.astype(np.float32)

    ss_res = np.sum((y_true - y_pred) ** 2, axis=0, dtype=np.float32)
    ss_tot = np.sum((y_true - np.mean(y_true, axis=0)) ** 2, axis=0, dtype=np.float32)

    return 1.0 - (ss_res / (ss_tot + 1e-8))


def load_masked_runs(nii_files, mask_img):
    """
    Load and mask all runs as float32
    """
    data = []
    for f in nii_files:
        img = nib.load(f)
        masked = apply_mask(img, mask_img).astype(np.float32)
        data.append(masked)
    return data


def load_design_matrices(events_files, n_regressors, permute=False, seed=None):
    """
    Load design matrices from CSV files.
    Each CSV should contain the design matrix for one run.
    """
    rng = np.random.default_rng(seed)
    designs = []

    for f in events_files:
        # Load CSV file
        df = pd.read_csv(f)
        
        # Convert to numpy array and select first n_regressors columns
        X = df.values.astype(np.float32)
        
        # If the CSV has more columns than needed, take only the first n_regressors
        if X.shape[1] > n_regressors:
            X = X[:, :n_regressors]
        elif X.shape[1] < n_regressors:
            raise ValueError(f"CSV {f} has {X.shape[1]} columns but {n_regressors} regressors requested")

        if permute:
            idx = rng.permutation(X.shape[0])
            X = X[idx]

        designs.append(X)

    return designs


def fit_glm(X, Y):
    """
    OLS using normal equation
    Float32
    """
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    XtX = X.T @ X
    XtY = X.T @ Y

    beta = np.linalg.pinv(XtX) @ XtY
    return beta.astype(np.float32)


# ============================================================
# Cross-validation core
# ============================================================

def run_cv(
    nii_files,
    events_files,
    mask_img,
    n_regressors,
    splits,
    random_state,
    permute=False,
):
    """
    Run cross-validated GLM and return mean R² map (float32)
    """

    masked_data = load_masked_runs(nii_files, mask_img)
    designs = load_design_matrices(
        events_files,
        n_regressors,
        permute=permute,
        seed=random_state,
    )

    coef_sum = None
    coef_sq_sum = None

    for split_id, (train_idx, test_idx) in enumerate(splits):

        # Concatenate train
        X_train = np.vstack([designs[i] for i in train_idx]).astype(np.float32)
        Y_train = np.vstack([masked_data[i] for i in train_idx]).astype(np.float32)

        beta = fit_glm(X_train, Y_train)

        # Concatenate test
        X_test = np.vstack([designs[i] for i in test_idx]).astype(np.float32)
        Y_test = np.vstack([masked_data[i] for i in test_idx]).astype(np.float32)

        Y_pred = X_test @ beta

        r2 = compute_r2(Y_test, Y_pred)

        if coef_sum is None:
            coef_sum = np.zeros_like(r2, dtype=np.float32)
            coef_sq_sum = np.zeros_like(r2, dtype=np.float32)

        coef_sum += r2
        coef_sq_sum += r2 ** 2

        del X_train, Y_train, X_test, Y_test, Y_pred, beta, r2
        gc.collect()

    n_splits = len(splits)

    mean_r2 = coef_sum / n_splits
    var_r2 = (coef_sq_sum / n_splits) - (mean_r2 ** 2)

    del coef_sum, coef_sq_sum
    gc.collect()

    return mean_r2.astype(np.float32), var_r2.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_regressors", type=int, required=True)
    parser.add_argument("--n_splits", type=int, default=50)
    parser.add_argument("--n_permutations", type=int, default=0)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--only_permutations", action="store_true")

    args = parser.parse_args()

    mask_img = nib.load(args.mask)

    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=0.5,
        random_state=args.random_state
    )

    splits = list(splitter.split(args.nii_files))

    # --------------------------------------------------------
    # REAL DATA
    # --------------------------------------------------------

    if not args.only_permutations:

        print("Running real CV...")
        mean_r2, var_r2 = run_cv(
            args.nii_files,
            args.events_files,
            mask_img,
            args.n_regressors,
            splits,
            args.random_state,
            permute=False
        )

        unmask(mean_r2, mask_img).to_filename(
            f"{args.output_prefix}_real_mean_r2.nii.gz"
        )

        unmask(var_r2, mask_img).to_filename(
            f"{args.output_prefix}_real_var_r2.nii.gz"
        )

        del mean_r2, var_r2
        gc.collect()

    # --------------------------------------------------------
    # PERMUTATIONS (Streaming, Memory Safe)
    # --------------------------------------------------------

    if args.n_permutations > 0:

        print(f"Running {args.n_permutations} permutations...")

        running_mean = None
        max_distribution = np.zeros(args.n_permutations, dtype=np.float32)

        for p in tqdm(range(args.n_permutations), desc="Permutations"):

            mean_r2, _ = run_cv(
                args.nii_files,
                args.events_files,
                mask_img,
                args.n_regressors,
                splits,
                args.random_state + p + 1,
                permute=True
            )

            if running_mean is None:
                running_mean = np.zeros_like(mean_r2, dtype=np.float32)

            running_mean += mean_r2
            max_distribution[p] = np.max(mean_r2)

            del mean_r2
            gc.collect()

        perm_mean = running_mean / args.n_permutations

        unmask(perm_mean, mask_img).to_filename(
            f"{args.output_prefix}_perm_mean_r2.nii.gz"
        )

        np.save(
            f"{args.output_prefix}_perm_max_distribution.npy",
            max_distribution
        )

        del perm_mean, running_mean
        gc.collect()

    print("Done.")


if __name__ == "__main__":
    main()
