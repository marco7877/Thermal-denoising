#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from nilearn.image import resample_to_img, concat_imgs, get_data
from nilearn.glm.first_level import FirstLevelModel
from sklearn.model_selection import ShuffleSplit
from tqdm import tqdm
import gc


# ============================================================
# Utilities
# ============================================================

def convert_to_percent_change(data, mask_img):
    """
    Convert time series to percent change: (x - mean)/mean * 100
    
    Parameters:
    -----------
    data : np.ndarray
        4D image data (x, y, z, time)
    mask_img : nibabel image
        Binary mask
    
    Returns:
    --------
    data_pc : np.ndarray
        Percent change data
    """
    mask_data = mask_img.get_fdata().astype(bool)
    
    # Reshape to voxels × time
    data_reshaped = data.reshape(-1, data.shape[-1])
    
    # Only compute for masked voxels to save memory
    masked_voxels = data_reshaped[mask_data.ravel(), :]
    
    # Compute mean for each voxel (across time)
    voxel_means = np.mean(masked_voxels, axis=1, keepdims=True)
    
    # Avoid division by zero
    voxel_means = np.where(voxel_means == 0, 1, voxel_means)
    
    # Convert to percent change
    masked_pc = (masked_voxels - voxel_means) / voxel_means * 100
    
    # Put back into full array
    data_pc_reshaped = np.zeros_like(data_reshaped)
    data_pc_reshaped[mask_data.ravel(), :] = masked_pc
    
    # Reshape back to 4D
    data_pc = data_pc_reshaped.reshape(data.shape)
    
    return data_pc


def regress_out_nuisance(data, design_matrix, n_task_regressors):
    """
    Regress out nuisance regressors from test time series data
    
    Parameters:
    -----------
    data : np.ndarray
        Time series data of shape (n_voxels, n_timepoints) - percent change
    design_matrix : pd.DataFrame
        Full design matrix (task + nuisance) for test runs
    n_task_regressors : int
        Number of task regressors (first columns)
    
    Returns:
    --------
    data_resid : np.ndarray
        Data with nuisance effects removed, same shape as input
    """
    # Get nuisance regressors (all columns after task regressors)
    nuisance_regressors = design_matrix.values[:, n_task_regressors:].astype(np.float32)
    
    # Add small regularization for numerical stability
    try:
        # Using pseudo-inverse
        pinv_nuisance = np.linalg.pinv(nuisance_regressors)
        # Projection matrix: H = X @ pinv(X)
        H_nuisance = nuisance_regressors @ pinv_nuisance
    except np.linalg.LinAlgError:
        # Fallback to regular inverse with regularization
        XtX = nuisance_regressors.T @ nuisance_regressors
        XtX_reg = XtX + 1e-8 * np.eye(XtX.shape[0])
        XtX_inv = np.linalg.inv(XtX_reg)
        H_nuisance = nuisance_regressors @ XtX_inv @ nuisance_regressors.T
    
    # Apply projection to remove nuisance effects
    # data shape: (n_voxels, n_timepoints)
    # H_nuisance shape: (n_timepoints, n_timepoints)
    # We want: data_cleaned = data - (data projected onto nuisance space)
    # Which is: data_cleaned = data - data @ H_nuisance
    data_resid = data - data @ H_nuisance
    
    return data_resid


def load_and_preprocess_runs(nii_files, mask_img):
    """
    Load all runs and convert to percent change
    Returns list of 4D arrays and the reference image
    """
    print("\nLoading and converting to percent change...")
    
    all_runs_data = []
    reference_img = None
    
    for i, f in enumerate(nii_files):
        print(f"  Processing run {i+1}/{len(nii_files)}: {f}")
        img = nib.load(f)
        
        # Use first run as reference
        if reference_img is None:
            reference_img = img
            print(f"    Using as reference")
        else:
            # Check if affine matches reference
            if not np.allclose(img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
                print(f"    Resampling to match reference space...")
                img = resample_to_img(img, reference_img, interpolation='continuous')
        
        # Get data and convert to percent change
        data = img.get_fdata().astype(np.float32)
        data_pc = convert_to_percent_change(data, mask_img)
        all_runs_data.append(data_pc)
        
        print(f"    Shape: {data_pc.shape}")
    
    return all_runs_data, reference_img


def load_design_matrices(events_files, n_task_regressors, permute=False, seed=None):
    """
    Load precomputed design matrices from CSV files
    These already have run-specific drift terms
    """
    rng = np.random.default_rng(seed)
    designs = []

    for i, f in enumerate(events_files):
        print(f"  Loading design matrix {i+1}/{len(events_files)}: {f}")
        
        # Load CSV file with apostrophe delimiter
        df = pd.read_csv(f, delimiter="'", quotechar=None, quoting=3, engine='python')
        
        # Remove any empty columns
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        
        if permute:
            # Permute the rows (timepoints) for null distribution
            df = df.sample(frac=1, random_state=rng).reset_index(drop=True)

        designs.append(df)

    return designs


def compute_r2_from_correlation(y_true, y_pred, mask):
    """
    Compute R² as squared Pearson correlation for masked voxels only
    """
    valid_voxels = mask.ravel() != 0
    
    y_true_reshaped = y_true.reshape(-1, y_true.shape[-1])
    y_pred_reshaped = y_pred.reshape(-1, y_pred.shape[-1])
    
    y_true_valid = y_true_reshaped[valid_voxels, :]
    y_pred_valid = y_pred_reshaped[valid_voxels, :]
    
    # Center the variables
    y_true_centered = y_true_valid - np.mean(y_true_valid, axis=1, keepdims=True)
    y_pred_centered = y_pred_valid - np.mean(y_pred_valid, axis=1, keepdims=True)
    
    # Compute Pearson correlation
    numerator = np.sum(y_true_centered * y_pred_centered, axis=1)
    denominator = np.sqrt(
        np.sum(y_true_centered**2, axis=1) * 
        np.sum(y_pred_centered**2, axis=1)
    )
    
    epsilon = 1e-8
    correlation = numerator / (denominator + epsilon)
    
    r2_values = correlation ** 2
    
    r2_full = np.zeros(y_true.shape[:-1], dtype=np.float32)
    r2_full[valid_voxels.reshape(y_true.shape[:-1])] = r2_values
    
    return r2_full


# ============================================================
# Cross-validation core
# ============================================================

def run_cv(
    nii_files,
    events_files,
    mask_img,
    n_task_regressors,
    tr,
    hrf_model,
    splits,
    random_state,
    permute=False,
):
    """
    Run cross-validated GLM with percent change data
    
    The logic:
    1. Training: Fit GLM on training runs using ALL regressors (task + nuisance)
    2. Extract task betas using contrast matrix (one beta per task regressor per voxel)
    3. Testing: 
       a. Regress out nuisance effects from test time series using test design matrix
       b. For each task, compute its contribution: beta_task * task_regressor_timecourse
       c. Sum all task contributions to get predicted time series
       d. Compute R² between cleaned test data and predicted time series
    """
    
    # Preprocess all runs to percent change
    all_runs_data, reference_img = load_and_preprocess_runs(nii_files, mask_img)
    
    # Load design matrices (already have run-specific drifts)
    print("\nLoading design matrices...")
    designs = load_design_matrices(
        events_files,
        n_task_regressors,
        permute=permute,
        seed=random_state,
    )
    
    # Get mask data
    mask_data = mask_img.get_fdata().astype(bool)
    mask_shape = mask_data.shape
    n_voxels_mask = np.sum(mask_data)
    print(f"\nMask has {n_voxels_mask} voxels")
    
    # Store results as 3D arrays (same as mask shape)
    r2_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sq_sum = np.zeros(mask_shape, dtype=np.float32)

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n{'='*60}")
        print(f"Processing split {split_id + 1}/{len(splits)}")
        print(f"{'='*60}")
        print(f"  Training runs: {train_idx}")
        print(f"  Testing runs: {test_idx}")
        
        # --- TRAINING ---
        # Get training data and design matrices
        train_data = [all_runs_data[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        
        # Concatenate training data along time dimension
        train_data_concat = np.concatenate(train_data, axis=-1)
        print(f"  Training data shape: {train_data_concat.shape}")
        
        # Concatenate design matrices
        train_design_matrix = pd.concat(train_designs, ignore_index=True)
        train_design_matrix = train_design_matrix.fillna(0)
        print(f"  Training design matrix shape: {train_design_matrix.shape}")
        
        # Create Nifti image for training data
        train_img = nib.Nifti1Image(train_data_concat, reference_img.affine)
        
        # Fit GLM
        print(f"  Fitting GLM on training data...")
        fmri_glm = FirstLevelModel(
            t_r=tr,
            mask_img=mask_img,
            standardize=False,
            signal_scaling=False,
            hrf_model=hrf_model,
            minimize_memory=True
        )
        
        fmri_glm = fmri_glm.fit(train_img, design_matrices=train_design_matrix)
        
        # Extract task betas (first n_task_regressors)
        all_regressors = train_design_matrix.columns.tolist()
        contrast_matrix = np.zeros((n_task_regressors, len(all_regressors)))
        for i in range(n_task_regressors):
            contrast_matrix[i, i] = 1
        
        print(f"  Extracting task betas...")
        betas_img = fmri_glm.compute_contrast(
            contrast_matrix,
            output_type='effect_size'
        )
        betas = get_data(betas_img).astype(np.float32)  # Shape: (x, y, z, n_task)
        print(f"  Betas shape: {betas.shape}")
        
        # --- TESTING ---
        # Get test data
        test_data_list = [all_runs_data[i] for i in test_idx]
        test_data = np.concatenate(test_data_list, axis=-1)  # Shape: (x, y, z, time)
        print(f"  Test data shape: {test_data.shape}")
        
        # Get test design matrices
        test_designs = [designs[i] for i in test_idx]
        test_design_matrix = pd.concat(test_designs, ignore_index=True)
        test_design_matrix = test_design_matrix.fillna(0)
        print(f"  Test design matrix shape: {test_design_matrix.shape}")
        
        # Reshape test data to (n_voxels, n_timepoints) for regression
        test_data_reshaped = test_data.reshape(-1, test_data.shape[-1])
        test_data_masked = test_data_reshaped[mask_data.ravel(), :]  # Shape: (n_voxels_mask, n_timepoints)
        print(f"  Test data masked shape: {test_data_masked.shape}")
        
        # CRITICAL STEP 1: Regress out nuisance effects from test data
        print(f"  Regressing out nuisance regressors from test data...")
        test_data_cleaned = regress_out_nuisance(
            test_data_masked,
            test_design_matrix,
            n_task_regressors
        )  # Shape: (n_voxels_mask, n_timepoints)
        print(f"  Cleaned test data shape: {test_data_cleaned.shape}")
        
        # CRITICAL STEP 2: Get task regressors for prediction
        # Each task regressor is a timecourse of that task's occurrence (convolved with HRF)
        test_task_regressors = test_design_matrix.values[:, :n_task_regressors].astype(np.float32)
        print(f"  Test task regressors shape: {test_task_regressors.shape}")  # (n_timepoints, n_task)
        
        # CRITICAL STEP 3: Reshape betas to (n_voxels_mask, n_task)
        betas_reshaped = betas.reshape(-1, n_task_regressors)
        betas_masked = betas_reshaped[mask_data.ravel(), :]  # Shape: (n_voxels_mask, n_task)
        print(f"  Betas masked shape: {betas_masked.shape}")
        
        # CRITICAL STEP 4: Predict by summing task contributions
        # For each task: beta_voxel * task_regressor_timecourse
        # Then sum across tasks
        print(f"  Predicting test time series (sum of task contributions)...")
        
        # Initialize prediction array
        predicted = np.zeros_like(test_data_cleaned)  # Shape: (n_voxels_mask, n_timepoints)
        
        # For each task, add its contribution
        for task_idx in range(n_task_regressors):
            task_beta = betas_masked[:, task_idx:task_idx+1]  # Shape: (n_voxels_mask, 1)
            task_regressor = test_task_regressors[:, task_idx:task_idx+1].T  # Shape: (1, n_timepoints)
            task_contribution = task_beta @ task_regressor  # Shape: (n_voxels_mask, n_timepoints)
            predicted += task_contribution
            print(f"    Task {task_idx+1}: contribution shape {task_contribution.shape}, "
                  f"beta range [{np.min(task_beta):.4f}, {np.max(task_beta):.4f}]")
        
        print(f"  Final predicted shape: {predicted.shape}")
        
        # Verify shapes match
        if test_data_cleaned.shape != predicted.shape:
            print(f"  WARNING: Shape mismatch!")
            print(f"    Cleaned test data: {test_data_cleaned.shape}")
            print(f"    Predicted: {predicted.shape}")
            # Take minimum shape
            min_voxels = min(test_data_cleaned.shape[0], predicted.shape[0])
            min_time = min(test_data_cleaned.shape[1], predicted.shape[1])
            test_data_cleaned = test_data_cleaned[:min_voxels, :min_time]
            predicted = predicted[:min_voxels, :min_time]
        
        # Compute R² between cleaned test data and predicted time series
        print(f"  Computing R²...")
        
        # Reconstruct 4D for R² computation
        test_data_4d = np.zeros(mask_shape + (test_data_cleaned.shape[1],), dtype=np.float32)
        predicted_4d = np.zeros(mask_shape + (predicted.shape[1],), dtype=np.float32)
        
        test_data_4d[mask_data, :] = test_data_cleaned
        predicted_4d[mask_data, :] = predicted
        
        r2_map = compute_r2_from_correlation(test_data_4d, predicted_4d, mask_data)
        
        # Calculate statistics
        r2_values = r2_map[mask_data]
        mean_r2 = np.mean(r2_values)
        std_r2 = np.std(r2_values)
        print(f"  R² statistics (within mask):")
        print(f"    Mean: {mean_r2:.6f}")
        print(f"    Std: {std_r2:.6f}")
        print(f"    Min: {np.min(r2_values):.6f}")
        print(f"    Max: {np.max(r2_values):.6f}")

        # Accumulate results (r2_map is already 3D)
        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2

        # Clean up
        del fmri_glm, betas, betas_reshaped, betas_masked, predicted, test_data, test_data_cleaned, r2_map
        gc.collect()

    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits  # This is 3D: (x, y, z)
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)

    print(f"\n{'='*60}")
    print(f"Final Results")
    print(f"{'='*60}")
    final_mean = np.mean(mean_r2[mask_data])
    final_std = np.std(mean_r2[mask_data])
    print(f"Mean R² across all splits (within mask): {final_mean:.6f}")
    print(f"Std R² across all splits: {final_std:.6f}")
    print(f"R² range: [{np.min(mean_r2[mask_data]):.6f}, {np.max(mean_r2[mask_data]):.6f}]")

    return mean_r2.astype(np.float32), var_r2.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Cross-validated GLM with percent change data")
    
    parser.add_argument("--nii_files", nargs="+", required=True,
                        help="List of NIfTI files (one per run)")
    parser.add_argument("--events_files", nargs="+", required=True,
                        help="List of design matrix CSV files (one per run)")
    parser.add_argument("--mask", required=True,
                        help="Mask NIfTI file")
    parser.add_argument("--n_task_regressors", type=int, required=True,
                        help="Number of task-related regressors (first columns in design matrix)")
    parser.add_argument("--tr", type=float, default=2.0,
                        help="Repetition time in seconds")
    parser.add_argument("--hrf_model", type=str, default="spm",
                        help="HRF model to use")
    parser.add_argument("--n_splits", type=int, default=50,
                        help="Number of cross-validation splits")
    parser.add_argument("--test_size", type=float, default=0.3,
                        help="Proportion of runs to use for testing")
    parser.add_argument("--n_permutations", type=int, default=0,
                        help="Number of permutations for null distribution")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random seed")
    parser.add_argument("--output_prefix", required=True,
                        help="Prefix for output files")
    parser.add_argument("--only_permutations", action="store_true",
                        help="Only run permutations (skip real data)")

    args = parser.parse_args()

    if len(args.nii_files) != len(args.events_files):
        raise ValueError(f"Number of NIfTI files ({len(args.nii_files)}) does not match number of event files ({len(args.events_files)})")

    print(f"\n{'='*60}")
    print(f"Cross-validated GLM with Percent Change Data")
    print(f"{'='*60}")
    print(f"Input files: {len(args.nii_files)} runs")
    print(f"Parameters:")
    print(f"  Task regressors: {args.n_task_regressors}")
    print(f"  TR: {args.tr}s")
    print(f"  HRF model: {args.hrf_model}")
    print(f"  Number of splits: {args.n_splits}")
    print(f"  Test size: {args.test_size*100:.0f}%")
    print(f"{'='*60}\n")

    # Load mask
    print("Loading mask...")
    mask_img = nib.load(args.mask)
    print(f"Mask shape: {mask_img.shape}")

    # Create cross-validation splits
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))

    # Run real data CV
    if not args.only_permutations:
        print("\n" + "="*60)
        print("RUNNING REAL DATA CV")
        print("="*60)
        
        mean_r2, var_r2 = run_cv(
            args.nii_files,
            args.events_files,
            mask_img,
            args.n_task_regressors,
            args.tr,
            args.hrf_model,
            splits,
            args.random_state,
            permute=False
        )

        # Save results - mean_r2 is 3D (x, y, z), which is exactly what unmask expects
        print("\nSaving results...")
        
        out_mean = f"{args.output_prefix}_real_mean_r2.nii.gz"
        # mean_r2 is already 3D with same shape as mask
        img_mean = nib.Nifti1Image(mean_r2, mask_img.affine)
        img_mean.to_filename(out_mean)
        print(f"  Saved: {out_mean}")

        out_var = f"{args.output_prefix}_real_var_r2.nii.gz"
        img_var = nib.Nifti1Image(var_r2, mask_img.affine)
        img_var.to_filename(out_var)
        print(f"  Saved: {out_var}")

        del mean_r2, var_r2
        gc.collect()

    # Run permutations
    if args.n_permutations > 0:
        print(f"\n{'='*60}")
        print(f"RUNNING {args.n_permutations} PERMUTATIONS")
        print(f"{'='*60}")

        running_mean = None
        max_distribution = np.zeros(args.n_permutations, dtype=np.float32)
        mask_data = mask_img.get_fdata().astype(bool)

        for p in tqdm(range(args.n_permutations), desc="Permutations"):
            mean_r2, _ = run_cv(
                args.nii_files,
                args.events_files,
                mask_img,
                args.n_task_regressors,
                args.tr,
                args.hrf_model,
                splits,
                args.random_state + p + 1,
                permute=True
            )

            if running_mean is None:
                running_mean = np.zeros_like(mean_r2, dtype=np.float32)

            running_mean += mean_r2
            max_distribution[p] = np.max(mean_r2[mask_data])

            del mean_r2
            gc.collect()

        # Save permutation results
        perm_mean = running_mean / args.n_permutations
        out_perm = f"{args.output_prefix}_perm_mean_r2.nii.gz"
        img_perm = nib.Nifti1Image(perm_mean, mask_img.affine)
        img_perm.to_filename(out_perm)
        print(f"  Saved: {out_perm}")

        out_dist = f"{args.output_prefix}_perm_max_distribution.npy"
        np.save(out_dist, max_distribution)
        print(f"  Saved: {out_dist}")

        # Calculate significance thresholds
        threshold_95 = np.percentile(max_distribution, 95)
        threshold_99 = np.percentile(max_distribution, 99)
        print(f"\nPermutation thresholds:")
        print(f"  95th: {threshold_95:.6f}")
        print(f"  99th: {threshold_99:.6f}")

        del perm_mean, running_mean
        gc.collect()

    print(f"\n{'='*60}")
    print("DONE!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
