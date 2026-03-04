#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from nilearn.image import resample_to_img, concat_imgs, get_data, load_img, clean_img, new_img_like
from nilearn.glm.first_level import FirstLevelModel
from sklearn.model_selection import ShuffleSplit
from tqdm import tqdm
import gc


# ============================================================
# Utilities
# ============================================================

def percent_change_img(img, mask_img):
    """
    Convert image to percent change: (x - mean)/mean * 100
    
    This is a simpler implementation that doesn't rely on internal nilearn functions
    """
    # Get data as numpy array
    data = get_data(img).astype(np.float64)  # Use float64 for precision
    mask_data = mask_img.get_fdata().astype(bool)
    
    # Get original shape
    original_shape = data.shape
    
    # Reshape to voxels × time
    n_voxels = np.prod(original_shape[:-1])
    n_timepoints = original_shape[-1]
    data_reshaped = data.reshape(n_voxels, n_timepoints)
    mask_flat = mask_data.ravel()
    
    print(f"    Debug - Original data range: [{np.min(data):.2f}, {np.max(data):.2f}]")
    print(f"    Debug - Original data mean: {np.mean(data):.2f}")
    
    # Only process masked voxels
    masked_data = data_reshaped[mask_flat, :]
    
    if masked_data.size > 0:
        # Compute mean for each masked voxel
        voxel_means = np.mean(masked_data, axis=1, keepdims=True)
        
        print(f"    Debug - Voxel means range: [{np.min(voxel_means):.2f}, {np.max(voxel_means):.2f}]")
        
        # Avoid division by zero
        voxel_means = np.where(voxel_means == 0, 1, voxel_means)
        
        # Convert to percent change
        masked_pc = (masked_data - voxel_means) / voxel_means * 100
        
        print(f"    Debug - Percent change range: [{np.min(masked_pc):.2f}, {np.max(masked_pc):.2f}]")
        
        # Put back into full array
        data_pc_reshaped = np.zeros_like(data_reshaped)
        data_pc_reshaped[mask_flat, :] = masked_pc
    else:
        data_pc_reshaped = data_reshaped
    
    # Reshape back to original 4D
    data_pc = data_pc_reshaped.reshape(original_shape).astype(np.float32)
    
    # Return as new image with same header/affine
    return new_img_like(img, data_pc)


def load_and_preprocess_runs(nii_files, mask_img):
    """
    Load all runs and convert to percent change FIRST
    Then they're ready for nuisance regression later
    """
    print("\nLoading runs and converting to percent change...")
    
    all_runs_imgs = []
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
        
        # CRITICAL STEP: Convert to percent change immediately
        print(f"    Converting to percent change...")
        img_pc = percent_change_img(img, mask_img)
        all_runs_imgs.append(img_pc)
        
        # Quick check of the data range
        data_sample = get_data(img_pc)[mask_img.get_fdata().astype(bool)]
        if len(data_sample) > 0:
            print(f"    Percent change range: [{np.min(data_sample):.2f}, {np.max(data_sample):.2f}]")
            print(f"    Percent change mean: {np.mean(data_sample):.2f}")
            print(f"    Percent change std: {np.std(data_sample):.2f}")
    
    return all_runs_imgs, reference_img


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
    Run cross-validated GLM with proper order:
    1. Convert ALL data to percent change FIRST
    2. Training: Fit GLM on training runs using ALL regressors
    3. Extract task betas
    4. Testing: 
       a. Regress out nuisance from percent-change test data
       b. Predict using task betas and task regressors
       c. Compute R²
    """
    
    # Step 1: Load all runs and convert to percent change FIRST
    all_runs_imgs, reference_img = load_and_preprocess_runs(nii_files, mask_img)
    
    # Load design matrices
    print("\nLoading design matrices...")
    designs = load_design_matrices(
        events_files,
        n_task_regressors,
        permute=permute,
        seed=random_state,
    )
    
    # Get mask data for later use
    mask_data = mask_img.get_fdata().astype(bool)
    mask_shape = mask_data.shape
    
    # Store results
    r2_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sq_sum = np.zeros(mask_shape, dtype=np.float32)

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n{'='*60}")
        print(f"Processing split {split_id + 1}/{len(splits)}")
        print(f"{'='*60}")
        print(f"  Training runs: {train_idx}")
        print(f"  Testing runs: {test_idx}")
        
        # --- TRAINING ---
        # Get training images (already in percent change) and design matrices
        train_imgs = [all_runs_imgs[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        
        # Concatenate training images
        print(f"  Concatenating training runs...")
        train_imgs_concat = concat_imgs(train_imgs)
        
        # Concatenate design matrices
        train_design_matrix = pd.concat(train_designs, ignore_index=True)
        train_design_matrix = train_design_matrix.fillna(0)
        print(f"  Training design matrix shape: {train_design_matrix.shape}")
        
        # Fit GLM on percent-change data
        print(f"  Fitting GLM on training data...")
        fmri_glm = FirstLevelModel(
            t_r=tr,
            mask_img=mask_img,
            standardize=False,  # Data already in percent change
            signal_scaling=False,  # Don't scale further
            hrf_model=hrf_model,
            minimize_memory=True
        )
        
        fmri_glm = fmri_glm.fit(train_imgs_concat, design_matrices=train_design_matrix)
        
        # Extract task betas
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
        print(f"  Betas range: [{np.min(betas):.4f}, {np.max(betas):.4f}]")
        
        # --- TESTING ---
        # Get test images (already in percent change) and design matrices
        test_imgs = [all_runs_imgs[i] for i in test_idx]
        test_designs = [designs[i] for i in test_idx]
        
        # Concatenate test images
        print(f"  Concatenating test runs...")
        test_imgs_concat = concat_imgs(test_imgs)
        
        # Get test design matrix
        test_design_matrix = pd.concat(test_designs, ignore_index=True)
        test_design_matrix = test_design_matrix.fillna(0)
        print(f"  Test design matrix shape: {test_design_matrix.shape}")
        
        # Step 2: Regress out nuisance regressors from percent-change test data
        print(f"  Regressing out nuisance regressors from test data...")
        
        # Extract nuisance regressors (all columns after task regressors)
        nuisance_regressors = test_design_matrix.values[:, n_task_regressors:].astype(np.float32)
        
        # Create confounds DataFrame
        confounds_df = pd.DataFrame(
            nuisance_regressors,
            columns=[f"nuisance_{i}" for i in range(nuisance_regressors.shape[1])]
        )
        
        # Use clean_img to remove nuisance effects from percent-change data
        test_imgs_cleaned = clean_img(
            test_imgs_concat,
            confounds=confounds_df,
            detrend=False,  # Drift terms already in confounds
            standardize=False,  # Keep percent change scale
            t_r=tr
        )
        print(f"  Cleaned test images shape: {test_imgs_cleaned.shape}")
        
        # Get the cleaned data
        test_data_cleaned = get_data(test_imgs_cleaned).astype(np.float32)
        print(f"  Cleaned test data range: [{np.min(test_data_cleaned):.2f}, {np.max(test_data_cleaned):.2f}]")
        print(f"  Cleaned test data mean: {np.mean(test_data_cleaned):.2f}")
        
        # Step 3: Get task regressors for prediction
        test_task_regressors = test_design_matrix.values[:, :n_task_regressors].astype(np.float32)
        print(f"  Test task regressors shape: {test_task_regressors.shape}")
        print(f"  Task regressors range: [{np.min(test_task_regressors):.4f}, {np.max(test_task_regressors):.4f}]")
        
        # Step 4: Reshape betas for prediction
        betas_reshaped = betas.reshape(-1, n_task_regressors)
        mask_flat = mask_data.ravel()
        betas_masked = betas_reshaped[mask_flat, :]
        print(f"  Betas masked shape: {betas_masked.shape}")
        print(f"  Betas masked range: [{np.min(betas_masked):.4f}, {np.max(betas_masked):.4f}]")
        
        # Step 5: Predict by summing task contributions
        print(f"  Predicting test time series...")
        
        # Reshape test data for comparison
        test_data_reshaped = test_data_cleaned.reshape(-1, test_data_cleaned.shape[-1])
        test_data_masked = test_data_reshaped[mask_flat, :]
        
        # Initialize prediction
        predicted = np.zeros_like(test_data_masked)
        
        # Sum task contributions
        for task_idx in range(n_task_regressors):
            task_beta = betas_masked[:, task_idx:task_idx+1]
            task_regressor = test_task_regressors[:, task_idx:task_idx+1].T
            task_contribution = task_beta @ task_regressor
            predicted += task_contribution
            print(f"    Task {task_idx+1} contribution range: [{np.min(task_contribution):.2f}, {np.max(task_contribution):.2f}]")
        
        print(f"  Predicted range: [{np.min(predicted):.2f}, {np.max(predicted):.2f}]")
        print(f"  Predicted mean: {np.mean(predicted):.2f}")
        print(f"  Cleaned test data range: [{np.min(test_data_masked):.2f}, {np.max(test_data_masked):.2f}]")
        print(f"  Cleaned test data mean: {np.mean(test_data_masked):.2f}")
        
        # Step 6: Compute R²
        print(f"  Computing R²...")
        
        # Reconstruct 4D
        test_data_4d = np.zeros(mask_shape + (test_data_masked.shape[1],), dtype=np.float32)
        predicted_4d = np.zeros(mask_shape + (predicted.shape[1],), dtype=np.float32)
        
        test_data_4d[mask_data, :] = test_data_masked
        predicted_4d[mask_data, :] = predicted
        
        r2_map = compute_r2_from_correlation(test_data_4d, predicted_4d, mask_data)
        
        # Statistics
        r2_values = r2_map[mask_data]
        mean_r2 = np.mean(r2_values)
        std_r2 = np.std(r2_values)
        print(f"  R² statistics (within mask):")
        print(f"    Mean: {mean_r2:.6f}")
        print(f"    Std: {std_r2:.6f}")
        print(f"    Min: {np.min(r2_values):.6f}")
        print(f"    Max: {np.max(r2_values):.6f}")

        # Accumulate
        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2

        # Clean up
        del fmri_glm, betas, betas_reshaped, betas_masked, predicted, test_data_cleaned, r2_map
        del test_imgs_cleaned, train_imgs_concat, test_imgs_concat
        gc.collect()

    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
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
    parser = argparse.ArgumentParser(description="Cross-validated GLM with percent change first")
    
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
    print(f"Cross-validated GLM - Percent Change First")
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

        # Save results
        print("\nSaving results...")
        
        out_mean = f"{args.output_prefix}_real_mean_r2.nii.gz"
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
