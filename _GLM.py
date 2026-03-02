#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from nilearn.image import resample_to_img, concat_imgs, get_data, load_img
from nilearn.glm.first_level import FirstLevelModel
from sklearn.model_selection import ShuffleSplit
from tqdm import tqdm
import gc
import warnings


# ============================================================
# Utilities
# ============================================================

def regress_out_nuisance(data, design_matrix, n_task_regressors):
    """
    Regress out nuisance regressors from time series data
    
    Parameters:
    -----------
    data : np.ndarray
        Time series data of shape (n_voxels, n_timepoints)
    design_matrix : pd.DataFrame
        Full design matrix (task + nuisance) - already includes constant
    n_task_regressors : int
        Number of task regressors (first columns)
    
    Returns:
    --------
    data_resid : np.ndarray
        Data with nuisance effects removed
    """
    # Get nuisance regressors (all columns after task regressors)
    # These already include constant term if present in original design
    nuisance_regressors = design_matrix.values[:, n_task_regressors:].astype(np.float32)
    
    # Compute nuisance projection matrix
    # H_nuisance = X_nuisance @ (X_nuisance.T @ X_nuisance)^{-1} @ X_nuisance.T
    try:
        # Using pseudo-inverse for numerical stability
        pinv_nuisance = np.linalg.pinv(nuisance_regressors)
        H_nuisance = nuisance_regressors @ pinv_nuisance
    except np.linalg.LinAlgError:
        # Fallback to regular inverse with regularization
        XtX = nuisance_regressors.T @ nuisance_regressors
        XtX_inv = np.linalg.inv(XtX + 1e-10 * np.eye(XtX.shape[0]))
        H_nuisance = nuisance_regressors @ XtX_inv @ nuisance_regressors.T
    
    # Apply projection to remove nuisance effects: (I - H_nuisance) @ data
    data_resid = data -  data @  H_nuisance
    
    return data_resid


def compute_r2_from_correlation(y_true, y_pred, mask):
    """
    Compute R² as squared Pearson correlation for masked voxels only
    Follows the formula from your working script:
    
    r = sum((x - mx)(y - my)) / sqrt(sum((x-mx)^2) * sum((y-my)^2))
    R² = r²
    """
    # Mask out zero voxels (assuming mask is binary)
    valid_voxels = mask.ravel() != 0
    
    # Reshape to voxels × time
    y_true_reshaped = y_true.reshape(-1, y_true.shape[-1])
    y_pred_reshaped = y_pred.reshape(-1, y_pred.shape[-1])
    
    # Only use valid voxels
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
    
    # Avoid division by zero
    epsilon = 1e-8
    correlation = numerator / (denominator + epsilon)
    
    # R² is squared correlation
    r2_values = correlation ** 2
    
    # Create full array with original shape
    r2_full = np.zeros(y_true.shape[:-1], dtype=np.float32)
    r2_full[valid_voxels.reshape(y_true.shape[:-1])] = r2_values
    
    return r2_full


def resample_to_reference(img, reference_img, interpolation='continuous'):
    """
    Resample image to match reference image's space
    """
    if not np.allclose(img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
        print(f"      Resampling to reference space...")
        img = resample_to_img(img, reference_img, interpolation=interpolation)
    return img


def load_and_resample_runs(nii_files, reference_img=None, tolerance=1e-3):
    """
    Load all runs and resample them to a common reference space
    If reference_img is None, use the first run as reference
    Returns list of resampled images and the reference image
    """
    images = []
    
    # Load first image to use as reference if none provided
    if reference_img is None:
        reference_img = nib.load(nii_files[0])
        print(f"  Using first run as reference: {nii_files[0]}")
    
    for i, f in enumerate(nii_files):
        print(f"  Loading and resampling run {i+1}/{len(nii_files)}: {f}")
        img = nib.load(f)
        
        # Resample to reference space
        img_resampled = resample_to_reference(img, reference_img)
        images.append(img_resampled)
        
    return images, reference_img


def load_masked_runs(nii_files, mask_img, reference_img=None, tolerance=1e-3):
    """
    Load and mask all runs, ensuring they're in the same space
    Returns list of masked data arrays
    """
    # First, ensure all runs are in the same space
    images, _ = load_and_resample_runs(nii_files, reference_img, tolerance)
    
    # Now apply mask to each resampled image
    data = []
    for i, img in enumerate(images):
        # Check if mask needs resampling to match image space
        if not np.allclose(img.affine, mask_img.affine, rtol=tolerance, atol=tolerance):
            print(f"    Resampling mask to match run {i+1} space...")
            mask_resampled = resample_to_img(mask_img, img, interpolation='nearest')
        else:
            mask_resampled = mask_img
            
        # Apply mask
        masked = apply_mask(img, mask_resampled).astype(np.float32)
        data.append(masked)
        
    return data, images[0]  # Return reference image for later use


def load_design_matrices(events_files, n_task_regressors, permute=False, seed=None):
    """
    Load design matrices from CSV files
    Each CSV should contain the full design matrix (task + nuisance)
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
            # Permute the rows (timepoints) of the design matrix
            df = df.sample(frac=1, random_state=rng).reset_index(drop=True)

        designs.append(df)

    return designs


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
    Run cross-validated GLM following your working approach:
    
    1. Training: Fit GLM on training runs using ALL regressors
    2. Extract task betas using contrast matrix
    3. Testing: 
       a. Regress out nuisance regressors from test time series
       b. Predict using ONLY task regressors from test design matrix
       c. Compute R² between nuisance-regressed test data and prediction
    """
    
    # Load brain data and get reference image
    print("\nLoading brain data...")
    brain_data, reference_img = load_masked_runs(nii_files, mask_img)
    
    # Load design matrices
    print("\nLoading design matrices...")
    designs = load_design_matrices(
        events_files,
        n_task_regressors,
        permute=permute,
        seed=random_state,
    )
    
    # Get mask data for later
    mask_data = mask_img.get_fdata().astype(bool)
    n_voxels = np.sum(mask_data)  # Number of voxels in mask
    
    # Store results
    r2_sum = None
    r2_sq_sum = None

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n  Processing split {split_id + 1}/{len(splits)}")
        print(f"    Training runs: {train_idx}")
        print(f"    Testing runs: {test_idx}")
        
        # --- TRAINING ---
        # Get training files and design matrices
        train_files = [nii_files[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        
        # Load and resample training images to reference space
        print(f"      Loading and resampling training runs...")
        train_imgs_list, _ = load_and_resample_runs(train_files, reference_img)
        
        # Concatenate training images
        print(f"      Concatenating training runs...")
        train_imgs = concat_imgs(train_imgs_list)
        
        # Concatenate design matrices
        train_design_matrix = pd.concat(train_designs, ignore_index=True)
        train_design_matrix = train_design_matrix.fillna(0)  # Replace NaN with 0
        
        # Fit GLM on training data
        print(f"      Fitting GLM on training data...")
        fmri_glm = FirstLevelModel(
            t_r=tr,
            mask_img=mask_img,
            standardize=False,
            signal_scaling=False,
            hrf_model=hrf_model,
            minimize_memory=True  # Set to True for memory efficiency
        )
        
        fmri_glm = fmri_glm.fit(train_imgs, design_matrices=train_design_matrix)
        
        # Create contrast matrix to extract task betas
        all_regressors = train_design_matrix.columns.tolist()
        
        # Create contrast matrix: one contrast per task regressor
        contrast_matrix = np.zeros((n_task_regressors, len(all_regressors)))
        for i in range(n_task_regressors):
            contrast_matrix[i, i] = 1
        
        # Extract task betas (effect sizes)
        print(f"      Extracting task betas...")
        betas_img = fmri_glm.compute_contrast(
            contrast_matrix,
            output_type='effect_size'
        )
        betas = get_data(betas_img).astype(np.float32)  # Shape: (x, y, z, n_task)
        
        # --- TESTING ---
        # Get test design matrices and concatenate
        test_designs = [designs[i] for i in test_idx]
        test_design_matrix = pd.concat(test_designs, ignore_index=True)
        test_design_matrix = test_design_matrix.fillna(0)
        
        # Get actual test data (already loaded and masked)
        test_data_list = [brain_data[i] for i in test_idx]
        # brain_data[i] is already masked and has shape (n_timepoints_i, n_voxels)
        # We need to concatenate along time dimension
        test_data = np.concatenate(test_data_list, axis=0)  # Shape: (n_timepoints_total, n_voxels)
        test_data = test_data.T  # Transpose to (n_voxels, n_timepoints_total)
        
        # CRITICAL STEP: Regress out nuisance effects from test data
        print(f"      Regressing out nuisance regressors from test data...")
        test_data_cleaned = regress_out_nuisance(
            test_data, 
            test_design_matrix, 
            n_task_regressors
        )  # Shape: (n_voxels, n_timepoints)
        
        # Take only task regressors for prediction
        test_task_regressors = test_design_matrix.values[:, :n_task_regressors].astype(np.float32)  # Shape: (n_timepoints, n_task)
        
        # Predict using tensor product
        print(f"      Predicting test time series...")
        # Reshape betas to (n_voxels, n_task)
        betas_reshaped = betas.reshape(-1, n_task_regressors)[mask_data.ravel(), :]  # Only keep masked voxels
        
        # Predicted time series: (n_voxels, n_timepoints)
        predicted = betas_reshaped @ test_task_regressors.T  # Shape: (n_voxels, n_timepoints)
        
        # Compute R² between cleaned test data and prediction
        print(f"      Computing R²...")
        
        # Both test_data_cleaned and predicted are already in shape (n_voxels, n_timepoints)
        # We need to reshape them back to 4D space for the R² computation
        # Create empty 4D arrays
        test_data_4d = np.zeros(mask_data.shape + (test_data_cleaned.shape[1],), dtype=np.float32)
        predicted_4d = np.zeros(mask_data.shape + (predicted.shape[1],), dtype=np.float32)
        
        # Fill in the masked voxels
        test_data_4d[mask_data, :] = test_data_cleaned
        predicted_4d[mask_data, :] = predicted
        
        r2_map = compute_r2_from_correlation(
            test_data_4d,
            predicted_4d,
            mask_data
        )
        
        mean_r2 = np.mean(r2_map[mask_data])
        print(f"      Mean R² (within mask): {mean_r2:.6f}")

        # Accumulate results
        if r2_sum is None:
            r2_sum = np.zeros_like(r2_map, dtype=np.float32)
            r2_sq_sum = np.zeros_like(r2_map, dtype=np.float32)

        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2

        # Clean up
        del fmri_glm, betas, betas_reshaped, predicted, test_data, test_data_cleaned, r2_map
        gc.collect()

    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)

    print(f"\n  Final mean R² across all splits (within mask): {np.mean(mean_r2[mask_data]):.6f}")

    return mean_r2.astype(np.float32), var_r2.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser(description="Cross-validated GLM following working approach")
    
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
                        help="HRF model to use (spm, glover, etc.)")
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

    # Verify input files match
    if len(args.nii_files) != len(args.events_files):
        raise ValueError(f"Number of NIfTI files ({len(args.nii_files)}) does not match number of event files ({len(args.events_files)})")

    print(f"\n{'='*60}")
    print(f"Cross-validated GLM (following working approach)")
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

    # Create cross-validation splits
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))

    # --------------------------------------------------------
    # REAL DATA
    # --------------------------------------------------------
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
        unmask(mean_r2, mask_img).to_filename(out_mean)
        print(f"  Saved: {out_mean}")
        
        mask_data = mask_img.get_fdata().astype(bool)
        print(f"    Mean R² (within mask): {np.mean(mean_r2[mask_data]):.6f}")
        print(f"    Max R²: {np.max(mean_r2[mask_data]):.6f}")

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
        unmask(perm_mean, mask_img).to_filename(out_perm)
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
