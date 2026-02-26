#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from nilearn.image import resample_to_img
from sklearn.model_selection import ShuffleSplit
from tqdm import tqdm
import gc
import warnings


# ============================================================
# Utilities
# ============================================================

def compute_r2(y_true, y_pred):
    """
    Proper cross-validated R² based on squared Pearson correlation
    R² = corr(y_true, y_pred)^2
    
    This ensures R² is between 0 and 1
    """
    y_true = y_true.astype(np.float32)
    y_pred = y_pred.astype(np.float32)
    
    # Center the variables
    y_true_centered = y_true - np.mean(y_true, axis=0, keepdims=True)
    y_pred_centered = y_pred - np.mean(y_pred, axis=0, keepdims=True)
    
    # Compute Pearson correlation
    numerator = np.sum(y_true_centered * y_pred_centered, axis=0)
    denominator = np.sqrt(np.sum(y_true_centered**2, axis=0) * np.sum(y_pred_centered**2, axis=0))
    
    # Avoid division by zero (where there's no variance in true or predicted)
    denominator = np.maximum(denominator, np.finfo(np.float32).eps)
    
    # Correlation coefficient
    corr = numerator / denominator
    
    # R² is squared correlation
    r2 = corr ** 2
    
    # Clip to handle any tiny numerical errors that might push it just above 1
    r2 = np.minimum(r2, 1.0)
    
    return r2.astype(np.float32)


def load_masked_runs(nii_files, mask_img, tolerance=1e-4):
    """
    Load and mask all runs as float32, resampling if needed
    Also return the number of timepoints per run for later use
    """
    data = []
    timepoints_per_run = []
    mask_affine = mask_img.affine
    
    for i, f in enumerate(nii_files):
        print(f"  Loading run {i+1}/{len(nii_files)}: {f}")
        img = nib.load(f)
        
        # Check if affine is different (with tolerance)
        if not np.allclose(img.affine, mask_affine, rtol=tolerance, atol=tolerance):
            print(f"    Affine mismatch detected. Resampling to match mask space...")
            print(f"    Max difference: {np.max(np.abs(img.affine - mask_affine))}")
            
            # Resample image to mask space
            img = resample_to_img(img, mask_img, interpolation='continuous')
            print(f"    Resampling complete")
        
        # Apply mask
        masked = apply_mask(img, mask_img).astype(np.float32)
        data.append(masked)
        timepoints_per_run.append(masked.shape[0])
        
    print(f"  Loaded {len(data)} runs with timepoints: {timepoints_per_run}")
    return data, timepoints_per_run


def load_design_matrices(events_files, total_regressors, permute=False, seed=None):
    """
    Load design matrices from CSV files with apostrophe delimiters.
    Each CSV should contain the design matrix for one run.
    Also return the number of timepoints per run
    """
    rng = np.random.default_rng(seed)
    designs = []
    timepoints_per_run = []

    for i, f in enumerate(events_files):
        print(f"  Loading design matrix {i+1}/{len(events_files)}: {f}")
        
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
        timepoints_per_run.append(X.shape[0])

    return designs, timepoints_per_run


def fit_glm(X, Y):
    """
    OLS using normal equation
    Float32
    
    Returns beta coefficients
    """
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    # Add small regularization for numerical stability
    n_features = X.shape[1]
    reg = np.eye(n_features) * 1e-6
    
    XtX = X.T @ X + reg
    XtY = X.T @ Y

    beta = np.linalg.solve(XtX, XtY)  # More stable than pinv
    return beta.astype(np.float32)


def predict(X, beta):
    """Make predictions from design matrix and beta coefficients"""
    return X @ beta


def residualize(X, Z):
    """
    Remove the effect of Z from X
    Returns X orthogonalized with respect to Z
    
    X: matrix to residualize (n_samples, n_features)
    Z: nuisance matrix (n_samples, n_nuisance)
    
    Returns: X_resid (n_samples, n_features)
    """
    # Add small regularization
    n_nuisance = Z.shape[1]
    reg = np.eye(n_nuisance) * 1e-6
    
    # Projection matrix for nuisance space
    ZtZ_inv = np.linalg.solve(Z.T @ Z + reg, np.eye(n_nuisance))
    P = Z @ ZtZ_inv @ Z.T
    
    # Residualize
    X_resid = X - P @ X
    
    return X_resid.astype(np.float32)


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
    
    The approach:
    1. In training: residualize task regressors w.r.t. nuisance regressors
    2. Fit model on residualized task regressors
    3. In testing: apply same residualization to test task regressors
    4. Predict and compute R² (squared correlation)
    
    This ensures we measure only the unique contribution of task regressors
    beyond what nuisance regressors can explain.
    """
    
    total_regressors = n_task_regressors + n_nuisance_regressors
    
    # Load data with timepoint information
    print("\nLoading brain data...")
    masked_data, brain_timepoints = load_masked_runs(nii_files, mask_img)
    
    print("\nLoading design matrices...")
    designs, design_timepoints = load_design_matrices(
        events_files,
        total_regressors,
        permute=permute,
        seed=random_state,
    )
    
    # Verify that timepoints match between brain data and design matrices
    print("\nVerifying timepoint compatibility...")
    for i, (bt, dt) in enumerate(zip(brain_timepoints, design_timepoints)):
        if bt != dt:
            raise ValueError(f"Run {i}: Brain data has {bt} timepoints but design matrix has {dt} timepoints")
    
    print(f"  All timepoints match: {brain_timepoints}")

    # Store results across splits
    r2_sum = None
    r2_sq_sum = None

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n  Processing split {split_id + 1}/{len(splits)}")
        print(f"    Training runs: {train_idx}")
        print(f"    Testing runs: {test_idx}")
        
        # --- TRAINING ---
        # Concatenate training data
        X_train_full = np.vstack([designs[i] for i in train_idx]).astype(np.float32)
        Y_train = np.vstack([masked_data[i] for i in train_idx]).astype(np.float32)
        
        # Split into task and nuisance regressors
        X_train_task = X_train_full[:, :n_task_regressors]
        X_train_nuisance = X_train_full[:, n_task_regressors:]
        
        # Remove nuisance effects from task regressors in training
        X_train_task_resid = residualize(X_train_task, X_train_nuisance)
        
        # Fit model using residualized task regressors
        beta_task = fit_glm(X_train_task_resid, Y_train)
        
        # --- TESTING ---
        # Process each test run separately to handle different lengths
        all_Y_test = []
        all_Y_pred = []
        
        for test_run_idx in test_idx:
            # Get this test run's data
            X_test_full = designs[test_run_idx].astype(np.float32)
            Y_test = masked_data[test_run_idx].astype(np.float32)
            
            # Split test data
            X_test_task = X_test_full[:, :n_task_regressors]
            X_test_nuisance = X_test_full[:, n_task_regressors:]
            
            # Apply same residualization to test task regressors
            # This ensures we're using the same transformation as training
            X_test_task_resid = residualize(X_test_task, X_test_nuisance)
            
            # Predict
            Y_pred = predict(X_test_task_resid, beta_task)
            
            # Store
            all_Y_test.append(Y_test)
            all_Y_pred.append(Y_pred)
        
        # Concatenate all test runs
        Y_test_all = np.vstack(all_Y_test)
        Y_pred_all = np.vstack(all_Y_pred)
        
        # Compute R² (squared correlation)
        r2 = compute_r2(Y_test_all, Y_pred_all)
        
        # Print statistics
        mean_r2 = np.mean(r2)
        print(f"    Mean R² across voxels: {mean_r2:.6f}")
        print(f"    R² range: [{np.min(r2):.6f}, {np.max(r2):.6f}]")
        print(f"    Voxels with R² > 0.1: {np.sum(r2 > 0.1)}/{r2.size}")

        # Accumulate for mean and variance
        if r2_sum is None:
            r2_sum = np.zeros_like(r2, dtype=np.float32)
            r2_sq_sum = np.zeros_like(r2, dtype=np.float32)

        r2_sum += r2
        r2_sq_sum += r2 ** 2

        # Clean up
        del X_train_full, X_train_task, X_train_nuisance, X_train_task_resid
        del X_test_task, X_test_nuisance, X_test_task_resid
        del Y_train, Y_test_all, Y_pred_all, beta_task, r2
        del all_Y_test, all_Y_pred
        gc.collect()

    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)

    # Final sanity check
    print(f"\n  Final mean R² across all splits and voxels: {np.mean(mean_r2):.6f}")
    print(f"  Final R² range: [{np.min(mean_r2):.6f}, {np.max(mean_r2):.6f}]")

    del r2_sum, r2_sq_sum
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

    # Verify that number of input files match
    if len(args.nii_files) != len(args.events_files):
        raise ValueError(f"Number of NIfTI files ({len(args.nii_files)}) does not match number of event files ({len(args.events_files)})")

    print(f"\n{'='*60}")
    print(f"Cross-validated GLM with Nuisance Regression")
    print(f"{'='*60}")
    print(f"Input files:")
    print(f"  NIfTI files: {len(args.nii_files)} runs")
    print(f"  Design files: {len(args.events_files)} runs")
    print(f"  Mask: {args.mask}")
    print(f"\nParameters:")
    print(f"  Task regressors: {args.n_task_regressors}")
    print(f"  Nuisance regressors: {args.n_nuisance_regressors}")
    print(f"  Total regressors: {args.n_task_regressors + args.n_nuisance_regressors}")
    print(f"  Number of splits: {args.n_splits}")
    print(f"  Test size: {args.test_size*100:.0f}%")
    print(f"  Train size: {(1-args.test_size)*100:.0f}%")
    print(f"  Permutations: {args.n_permutations}")
    print(f"  Random seed: {args.random_state}")
    print(f"  Output prefix: {args.output_prefix}")
    print(f"{'='*60}\n")

    # Load mask
    print("Loading mask...")
    mask_img = nib.load(args.mask)
    print(f"  Mask shape: {mask_img.shape}")

    # Create cross-validation splits with specified test_size
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))
    
    print(f"\nCross-validation splits created: {len(splits)}")

    # --------------------------------------------------------
    # REAL DATA
    # --------------------------------------------------------
    if not args.only_permutations:
        print("\n" + "="*60)
        print("RUNNING REAL DATA CV WITH NUISANCE REGRESSION")
        print("="*60)
        
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
        print("\nSaving results...")
        
        out_mean = f"{args.output_prefix}_real_mean_r2.nii.gz"
        unmask(mean_r2, mask_img).to_filename(out_mean)
        print(f"  Saved: {out_mean}")
        print(f"    Mean R² across all voxels: {np.mean(mean_r2):.6f}")
        print(f"    Std R² across voxels: {np.std(mean_r2):.6f}")
        print(f"    Max R²: {np.max(mean_r2):.6f}")
        print(f"    Min R²: {np.min(mean_r2):.6f}")

        out_var = f"{args.output_prefix}_real_var_r2.nii.gz"
        unmask(var_r2, mask_img).to_filename(out_var)
        print(f"  Saved: {out_var}")

        del mean_r2, var_r2
        gc.collect()

    # --------------------------------------------------------
    # PERMUTATIONS
    # --------------------------------------------------------
    if args.n_permutations > 0:
        print(f"\n{'='*60}")
        print(f"RUNNING {args.n_permutations} PERMUTATIONS")
        print(f"{'='*60}")

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
                args.random_state + p + 1,
                permute=True
            )

            if running_mean is None:
                running_mean = np.zeros_like(mean_r2, dtype=np.float32)

            running_mean += mean_r2
            max_distribution[p] = np.max(mean_r2)

            del mean_r2
            gc.collect()

        # Save permutation results
        print("\nSaving permutation results...")
        
        perm_mean = running_mean / args.n_permutations
        out_perm = f"{args.output_prefix}_perm_mean_r2.nii.gz"
        unmask(perm_mean, mask_img).to_filename(out_perm)
        print(f"  Saved: {out_perm}")
        print(f"    Mean R² across permutations: {np.mean(perm_mean):.6f}")

        out_dist = f"{args.output_prefix}_perm_max_distribution.npy"
        np.save(out_dist, max_distribution)
        print(f"  Saved: {out_dist}")

        # Calculate and print significance threshold
        threshold_95 = np.percentile(max_distribution, 95)
        threshold_99 = np.percentile(max_distribution, 99)
        print(f"\nPermutation-based significance thresholds:")
        print(f"  95th percentile: {threshold_95:.6f}")
        print(f"  99th percentile: {threshold_99:.6f}")
        print(f"  Max value in permutations: {np.max(max_distribution):.6f}")

        del perm_mean, running_mean
        gc.collect()

    print(f"\n{'='*60}")
    print("DONE!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
