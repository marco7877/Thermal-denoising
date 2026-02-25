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


def load_design_matrices(events_files, total_regressors, permute=False, seed=None):
    """
    Load design matrices from CSV files with apostrophe delimiters.
    Each CSV should contain the design matrix for one run.
    """
    rng = np.random.default_rng(seed)
    designs = []

    for f in events_files:
        # Load CSV file with apostrophe delimiter
        df = pd.read_csv(f, delimiter="'", quotechar=None, quoting=3, engine='python')
        
        # Remove any empty columns that might result from the delimiter parsing
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        
        # Convert to numpy array and ensure numeric type
        X = df.values.astype(np.float32)
        
        # Check if we have enough columns
        if X.shape[1] > total_regressors:
            X = X[:, :total_regressors]
        elif X.shape[1] < total_regressors:
            raise ValueError(f"CSV {f} has {X.shape[1]} columns but {total_regressors} regressors requested")

        if permute:
            # For permutations, we permute the rows (timepoints) of the design matrix
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


def orthogonalize_regressors(X_task, X_nuisance):
    """
    Orthogonalize task regressors with respect to nuisance regressors
    Returns task regressors with nuisance variance removed
    """
    # Projection matrix for nuisance space
    X_nuisance = X_nuisance.astype(np.float32)
    X_task = X_task.astype(np.float32)
    
    # Compute projection matrix: H = X_nuisance @ pinv(X_nuisance.T @ X_nuisance) @ X_nuisance.T
    # But we can do it more efficiently:
    nuisance_params = np.linalg.pinv(X_nuisance.T @ X_nuisance) @ X_nuisance.T @ X_task
    X_task_predicted = X_nuisance @ nuisance_params
    
    # Residuals are task regressors with nuisance variance removed
    X_task_orth = X_task - X_task_predicted
    
    return X_task_orth.astype(np.float32)


# ============================================================
# Cross-validation core with nuisance regression
# ============================================================

def run_cv(
    nii_files,
    events_files,
    mask_img,
    n_task_regressors,
    n_nuisance_regressors,
    splits,
    random_state,
    permute=False,
):
    """
    Run cross-validated GLM with nuisance regression and return mean R² map (float32)
    
    Parameters:
    -----------
    n_task_regressors : int
        Number of task-related regressors (first columns in design matrix)
    n_nuisance_regressors : int
        Number of nuisance regressors (drifts, constant, etc.) - remaining columns
    """
    
    total_regressors = n_task_regressors + n_nuisance_regressors
    
    # Load data
    masked_data = load_masked_runs(nii_files, mask_img)
    designs = load_design_matrices(
        events_files,
        total_regressors,
        permute=permute,
        seed=random_state,
    )

    coef_sum = None
    coef_sq_sum = None

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        # --- TRAINING ---
        # Concatenate training data
        X_train_full = np.vstack([designs[i] for i in train_idx]).astype(np.float32)
        Y_train = np.vstack([masked_data[i] for i in train_idx]).astype(np.float32)
        
        # Split into task and nuisance regressors
        X_train_task = X_train_full[:, :n_task_regressors]
        X_train_nuisance = X_train_full[:, n_task_regressors:]
        
        # Orthogonalize task regressors with respect to nuisance in training set
        X_train_task_orth = orthogonalize_regressors(X_train_task, X_train_nuisance)
        
        # Fit model using orthogonalized task regressors
        beta_task = fit_glm(X_train_task_orth, Y_train)
        
        # --- TESTING ---
        # Concatenate test data
        X_test_full = np.vstack([designs[i] for i in test_idx]).astype(np.float32)
        Y_test = np.vstack([masked_data[i] for i in test_idx]).astype(np.float32)
        
        # Split test data
        X_test_task = X_test_full[:, :n_task_regressors]
        X_test_nuisance = X_test_full[:, n_task_regressors:]
        
        # Orthogonalize test task regressors using the same nuisance space (from training)
        # This ensures we don't use test data to estimate nuisance parameters
        nuisance_params = np.linalg.pinv(X_train_nuisance.T @ X_train_nuisance) @ X_train_nuisance.T @ X_test_task
        X_test_task_predicted = X_test_nuisance @ nuisance_params
        X_test_task_orth = X_test_task - X_test_task_predicted
        
        # Predict using task model
        Y_pred = X_test_task_orth @ beta_task
        
        # Compute R² - this now reflects variance explained by task regressors only
        r2 = compute_r2(Y_test, Y_pred)

        # Accumulate for mean and variance
        if coef_sum is None:
            coef_sum = np.zeros_like(r2, dtype=np.float32)
            coef_sq_sum = np.zeros_like(r2, dtype=np.float32)

        coef_sum += r2
        coef_sq_sum += r2 ** 2

        # Clean up
        del X_train_full, X_train_task, X_train_nuisance, X_train_task_orth
        del X_test_full, X_test_task, X_test_nuisance, X_test_task_orth
        del Y_train, Y_test, Y_pred, beta_task, r2, nuisance_params
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

    parser = argparse.ArgumentParser(description="Cross-validated GLM with nuisance regression")
    
    parser.add_argument("--nii_files", nargs="+", required=True,
                        help="List of NIfTI files (one per run)")
    parser.add_argument("--events_files", nargs="+", required=True,
                        help="List of design matrix CSV files (one per run)")
    parser.add_argument("--mask", required=True,
                        help="Mask NIfTI file")
    parser.add_argument("--n_task_regressors", type=int, required=True,
                        help="Number of task-related regressors (first columns in design matrix)")
    parser.add_argument("--n_nuisance_regressors", type=int, required=True,
                        help="Number of nuisance regressors (drifts, constant, etc.)")
    parser.add_argument("--n_splits", type=int, default=50,
                        help="Number of cross-validation splits")
    parser.add_argument("--test_size", type=float, default=0.5,
                        help="Proportion of runs to use for testing (default: 0.5)")
    parser.add_argument("--n_permutations", type=int, default=0,
                        help="Number of permutations for null distribution")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--output_prefix", required=True,
                        help="Prefix for output files")
    parser.add_argument("--only_permutations", action="store_true",
                        help="Only run permutations (skip real data)")

    args = parser.parse_args()

    # Load mask
    mask_img = nib.load(args.mask)

    # Create cross-validation splits with specified test_size
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))
    
    print(f"\nCross-validation setup:")
    print(f"  Number of splits: {args.n_splits}")
    print(f"  Test size: {args.test_size*100:.0f}%")
    print(f"  Train size: {(1-args.test_size)*100:.0f}%")

    # --------------------------------------------------------
    # REAL DATA
    # --------------------------------------------------------
    if not args.only_permutations:
        print("\n" + "="*50)
        print("Running real data CV with nuisance regression...")
        print(f"Task regressors: {args.n_task_regressors}")
        print(f"Nuisance regressors: {args.n_nuisance_regressors}")
        print("="*50)
        
        mean_r2, var_r2 = run_cv(
            args.nii_files,
            args.events_files,
            mask_img,
            args.n_task_regressors,
            args.n_nuisance_regressors,
            splits,
            args.random_state,
            permute=False
        )

        # Save results
        unmask(mean_r2, mask_img).to_filename(
            f"{args.output_prefix}_real_mean_r2.nii.gz"
        )
        print(f"Saved: {args.output_prefix}_real_mean_r2.nii.gz")

        unmask(var_r2, mask_img).to_filename(
            f"{args.output_prefix}_real_var_r2.nii.gz"
        )
        print(f"Saved: {args.output_prefix}_real_var_r2.nii.gz")

        del mean_r2, var_r2
        gc.collect()

    # --------------------------------------------------------
    # PERMUTATIONS
    # --------------------------------------------------------
    if args.n_permutations > 0:
        print(f"\nRunning {args.n_permutations} permutations...")
        print("="*50)

        running_mean = None
        max_distribution = np.zeros(args.n_permutations, dtype=np.float32)

        for p in tqdm(range(args.n_permutations), desc="Permutations"):
            mean_r2, _ = run_cv(
                args.nii_files,
                args.events_files,
                mask_img,
                args.n_task_regressors,
                args.n_nuisance_regressors,
                splits,
                args.random_state + p + 1,  # Different seed for each permutation
                permute=True
            )

            if running_mean is None:
                running_mean = np.zeros_like(mean_r2, dtype=np.float32)

            running_mean += mean_r2
            max_distribution[p] = np.max(mean_r2)

            del mean_r2
            gc.collect()

        # Save permutation results
        perm_mean = running_mean / args.n_permutations
        unmask(perm_mean, mask_img).to_filename(
            f"{args.output_prefix}_perm_mean_r2.nii.gz"
        )
        print(f"Saved: {args.output_prefix}_perm_mean_r2.nii.gz")

        np.save(
            f"{args.output_prefix}_perm_max_distribution.npy",
            max_distribution
        )
        print(f"Saved: {args.output_prefix}_perm_max_distribution.npy")

        # Calculate and print significance threshold
        threshold_95 = np.percentile(max_distribution, 95)
        threshold_99 = np.percentile(max_distribution, 99)
        print(f"\nPermutation-based significance thresholds:")
        print(f"95th percentile: {threshold_95:.4f}")
        print(f"99th percentile: {threshold_99:.4f}")

        del perm_mean, running_mean
        gc.collect()

    print("\nDone!")


if __name__ == "__main__":
    main()
