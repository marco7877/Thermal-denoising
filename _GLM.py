#!/usr/bin/env python

import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from nilearn.masking import apply_mask, unmask
from nilearn.image import resample_to_img, concat_imgs, get_data, load_img, clean_img, new_img_like
from nilearn.glm.first_level import FirstLevelModel
from sklearn.model_selection import ShuffleSplit
from sklearn.metrics import r2_score
from scipy import stats
from tqdm import tqdm
import gc
import warnings
warnings.filterwarnings('ignore')


# ============================================================
# Utilities
# ============================================================

def percent_change_scaling(data, mask_img):
    """
    Apply percent change scaling to match GLM's signal_scaling='psc'
    This is (x - mean)/mean * 100
    """
    mask_data = mask_img.get_fdata().astype(bool)
    
    # Reshape to voxels × time
    original_shape = data.shape
    data_reshaped = data.reshape(-1, original_shape[-1])
    mask_flat = mask_data.ravel()
    
    # Compute mean for each voxel
    voxel_means = np.mean(data_reshaped[mask_flat, :], axis=1, keepdims=True)
    
    # Avoid division by zero
    voxel_means = np.where(np.abs(voxel_means) < 1e-6, 1, voxel_means)
    
    # Apply percent change scaling to masked voxels
    data_scaled_reshaped = np.zeros_like(data_reshaped)
    data_scaled_reshaped[mask_flat, :] = (data_reshaped[mask_flat, :] - voxel_means) / voxel_means * 100
    
    return data_scaled_reshaped.reshape(original_shape)


def zscore_design_matrix(design_matrix):
    """
    Manually z-score the design matrix to match GLM's standardize=True
    """
    design_scaled = design_matrix.copy()
    for col in design_scaled.columns:
        mean = design_scaled[col].mean()
        std = design_scaled[col].std()
        if std > 0:
            design_scaled[col] = (design_scaled[col] - mean) / std
        else:
            design_scaled[col] = 0  # Constant column
    return design_scaled


def create_block_diagonal_design(designs, n_task_regressors):
    """
    Create a proper block-diagonal design matrix where each run has its own
    nuisance regressors (drifts, constant) that are independent across runs.
    
    This is CRITICAL: polynomials must be run-specific, not shared across runs.
    
    Parameters:
    -----------
    designs : list of pd.DataFrame
        List of design matrices, one per run
    n_task_regressors : int
        Number of task regressors (assumed to be the same across runs)
    
    Returns:
    --------
    block_design : pd.DataFrame
        Block-diagonal design matrix with run-specific nuisance regressors
    """
    all_blocks = []
    
    for run_idx, design in enumerate(designs):
        # Split into task and nuisance parts
        task_part = design.iloc[:, :n_task_regressors].copy()
        nuisance_part = design.iloc[:, n_task_regressors:].copy()
        
        # Rename nuisance regressors to be run-specific
        # This ensures they are treated independently across runs
        rename_dict = {}
        for col in nuisance_part.columns:
            rename_dict[col] = f"{col}_run{run_idx}"
        nuisance_part.rename(columns=rename_dict, inplace=True)
        
        # Combine task and renamed nuisance for this run
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        all_blocks.append(run_design)
    
    # Concatenate vertically (stack runs)
    block_design = pd.concat(all_blocks, axis=0, ignore_index=True)
    block_design = block_design.fillna(0)
    
    return block_design


def load_and_preprocess_runs(nii_files, mask_img):
    """
    Load all runs and return with their number of timepoints
    """
    print("\n" + "="*60)
    print("STEP 1: Loading runs")
    print("="*60)
    
    all_runs_imgs = []
    run_timepoints = []
    reference_img = None
    
    for i, f in enumerate(nii_files):
        print(f"\n  Processing run {i+1}/{len(nii_files)}: {f}")
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
        
        # Get number of timepoints
        n_timepoints = img.shape[-1]
        run_timepoints.append(n_timepoints)
        all_runs_imgs.append(img)
        
        print(f"    Shape: {img.shape}, Timepoints: {n_timepoints}")
    
    # Print summary
    print(f"\n  Summary of run timepoints:")
    for i, tp in enumerate(run_timepoints):
        print(f"    Run {i}: {tp} timepoints")
    print(f"    Total timepoints across all runs: {sum(run_timepoints)}")
    
    return all_runs_imgs, reference_img, run_timepoints


def load_design_matrices(events_files, run_timepoints, n_task_regressors, permute=False, seed=None):
    """
    Load precomputed design matrices and verify they match run timepoints
    """
    rng = np.random.default_rng(seed)
    designs = []

    print("\n" + "="*60)
    print("STEP 2: Loading design matrices")
    print("="*60)

    for i, f in enumerate(events_files):
        print(f"\n  Loading design matrix {i+1}/{len(events_files)}: {f}")
        
        # Load CSV file
        df = pd.read_csv(f)
        
        # Remove any empty columns
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        
        print(f"    Design matrix has {len(df)} rows")
        print(f"    fMRI run {i} has {run_timepoints[i]} timepoints")
        
        # Check if number of rows matches run timepoints
        if len(df) != run_timepoints[i]:
            print(f"    WARNING: Mismatch in timepoints for run {i}!")
            print(f"      Design matrix: {len(df)} rows")
            print(f"      fMRI run: {run_timepoints[i]} timepoints")
            raise ValueError(f"Timepoint mismatch for run {i}: fMRI={run_timepoints[i]}, Design={len(df)}")
        
        if permute:
            df = df.sample(frac=1, random_state=rng).reset_index(drop=True)

        designs.append(df)
        
        # Print column structure for first file
        if i == 0:
            print(f"    Columns: {df.columns.tolist()}")
            print(f"    Task regressors (first {n_task_regressors}): {df.columns[:n_task_regressors].tolist()}")
            print(f"    Nuisance regressors: {df.columns[n_task_regressors:].tolist()}")

    # Print summary
    print(f"\n  Summary of design matrix rows:")
    for i, design in enumerate(designs):
        print(f"    Run {i}: {len(design)} rows")
    print(f"    Total rows across all designs: {sum(len(d) for d in designs)}")

    return designs


def compute_r2_multiple_methods(y_true, y_pred, mask):
    """
    Compute R² using multiple methods for verification
    """
    valid_voxels = mask.ravel() != 0
    
    y_true_reshaped = y_true.reshape(-1, y_true.shape[-1])
    y_pred_reshaped = y_pred.reshape(-1, y_pred.shape[-1])
    
    y_true_valid = y_true_reshaped[valid_voxels, :]
    y_pred_valid = y_pred_reshaped[valid_voxels, :]
    
    n_voxels = y_true_valid.shape[0]
    n_timepoints = y_true_valid.shape[1]
    
    print(f"    R² computation: {n_voxels} voxels, {n_timepoints} timepoints")
    
    # Initialize arrays for each method
    r2_corr = np.zeros(n_voxels)
    r2_1minus = np.zeros(n_voxels)
    r2_sklearn = np.zeros(n_voxels)
    
    # Additional diagnostics
    mean_correlation = 0
    neg_count_1minus = 0
    neg_count_sklearn = 0
    large_diff_count = 0
    
    for v in range(n_voxels):
        y_t = y_true_valid[v, :]
        y_p = y_pred_valid[v, :]
        
        # Method 1: Squared Pearson correlation
        corr = np.corrcoef(y_t, y_p)[0, 1]
        if not np.isnan(corr):
            r2_corr[v] = corr ** 2
            mean_correlation += abs(corr)
        else:
            r2_corr[v] = 0
        
        # Method 2: 1 - (SS_residual / SS_total)
        ss_res = np.sum((y_t - y_p) ** 2)
        ss_tot = np.sum((y_t - np.mean(y_t)) ** 2)
        r2_1minus[v] = 1 - (ss_res / (ss_tot + 1e-8))
        if r2_1minus[v] < 0:
            neg_count_1minus += 1
        
        # Method 3: sklearn's r2_score
        r2_sklearn[v] = r2_score(y_t, y_p)
        if r2_sklearn[v] < 0:
            neg_count_sklearn += 1
        
        # Check for large discrepancies between methods
        if abs(r2_corr[v] - r2_1minus[v]) > 0.01:
            large_diff_count += 1
    
    mean_correlation /= n_voxels
    
    # Create full 3D maps
    r2_corr_full = np.zeros(y_true.shape[:-1], dtype=np.float32)
    r2_1minus_full = np.zeros(y_true.shape[:-1], dtype=np.float32)
    r2_sklearn_full = np.zeros(y_true.shape[:-1], dtype=np.float32)
    
    r2_corr_full[valid_voxels.reshape(y_true.shape[:-1])] = r2_corr
    r2_1minus_full[valid_voxels.reshape(y_true.shape[:-1])] = r2_1minus
    r2_sklearn_full[valid_voxels.reshape(y_true.shape[:-1])] = r2_sklearn
    
    # Print comparison of methods
    print(f"\n    {'='*50}")
    print(f"    R² METHOD COMPARISON")
    print(f"    {'='*50}")
    print(f"    Method 1 (Squared Correlation):")
    print(f"      Mean: {np.mean(r2_corr):.6f}")
    print(f"      Std:  {np.std(r2_corr):.6f}")
    print(f"      Min:  {np.min(r2_corr):.6f}")
    print(f"      Max:  {np.max(r2_corr):.6f}")
    
    print(f"\n    Method 2 (1 - SS_res/SS_tot):")
    print(f"      Mean: {np.mean(r2_1minus):.6f}")
    print(f"      Std:  {np.std(r2_1minus):.6f}")
    print(f"      Min:  {np.min(r2_1minus):.6f}")
    print(f"      Max:  {np.max(r2_1minus):.6f}")
    print(f"      Negative R² voxels: {neg_count_1minus}/{n_voxels} ({100*neg_count_1minus/n_voxels:.2f}%)")
    
    print(f"\n    Method 3 (sklearn r2_score):")
    print(f"      Mean: {np.mean(r2_sklearn):.6f}")
    print(f"      Std:  {np.std(r2_sklearn):.6f}")
    print(f"      Min:  {np.min(r2_sklearn):.6f}")
    print(f"      Max:  {np.max(r2_sklearn):.6f}")
    print(f"      Negative R² voxels: {neg_count_sklearn}/{n_voxels} ({100*neg_count_sklearn/n_voxels:.2f}%)")
    
    print(f"\n    Cross-method comparisons:")
    print(f"      Mean absolute correlation: {mean_correlation:.6f}")
    print(f"      Correlation between Method 1 and 2: {np.corrcoef(r2_corr, r2_1minus)[0,1]:.6f}")
    print(f"      Correlation between Method 1 and 3: {np.corrcoef(r2_corr, r2_sklearn)[0,1]:.6f}")
    print(f"      Voxels with >0.01 difference: {large_diff_count}/{n_voxels} ({100*large_diff_count/n_voxels:.2f}%)")
    print(f"    {'='*50}")
    
    return r2_corr_full, r2_1minus_full, r2_sklearn_full


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
    Run cross-validated GLM with proper block-diagonal design matrices
    """
    
    # Step 1: Load all runs and get timepoints
    all_runs_imgs, reference_img, run_timepoints = load_and_preprocess_runs(nii_files, mask_img)
    
    # Load design matrices with timepoint alignment
    designs = load_design_matrices(
        events_files,
        run_timepoints,
        n_task_regressors,
        permute=permute,
        seed=random_state,
    )
    
    # Get mask data
    mask_data = mask_img.get_fdata().astype(bool)
    mask_shape = mask_data.shape
    print(f"\nMask has {np.sum(mask_data)} voxels")
    
    # Store results for each method
    r2_corr_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_1minus_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sklearn_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_corr_sq_sum = np.zeros(mask_shape, dtype=np.float32)

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n{'='*60}")
        print(f"STEP 3: Processing split {split_id + 1}/{len(splits)}")
        print(f"{'='*60}")
        print(f"  Training runs: {train_idx}")
        print(f"  Testing runs: {test_idx}")
        
        # --- TRAINING ---
        print(f"\n  --- Training Phase ---")
        train_imgs = [all_runs_imgs[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        
        # Calculate total timepoints for training
        train_total_tp = sum(run_timepoints[i] for i in train_idx)
        train_design_total_tp = sum(len(designs[i]) for i in train_idx)
        print(f"  Training total timepoints: fMRI={train_total_tp}, Design={train_design_total_tp}")
        
        if train_total_tp != train_design_total_tp:
            raise ValueError(f"Training timepoint mismatch! fMRI={train_total_tp}, Design={train_design_total_tp}")
        
        print(f"  Concatenating {len(train_imgs)} training runs...")
        train_imgs_concat = concat_imgs(train_imgs)
        train_data = get_data(train_imgs_concat)
        print(f"    Training data shape: {train_data.shape}")
        print(f"    Training data range: [{np.min(train_data):.2f}, {np.max(train_data):.2f}]")
        
        # Apply percent change scaling to training data
        print(f"  Applying percent change scaling to training data...")
        train_data_scaled = percent_change_scaling(train_data, mask_img)
        train_img_scaled = new_img_like(train_imgs_concat, train_data_scaled)
        print(f"    Scaled training data range: [{np.min(train_data_scaled):.2f}, {np.max(train_data_scaled):.2f}]")
        
        # CRITICAL: Create block-diagonal design matrix with run-specific nuisance regressors
        print(f"  Creating block-diagonal design matrix with run-specific nuisance regressors...")
        train_design_block = create_block_diagonal_design(train_designs, n_task_regressors)
        print(f"    Block design shape: {train_design_block.shape}")
        print(f"    Columns: {train_design_block.columns.tolist()}")
        
        # Z-score the design matrix
        train_design_zscored = zscore_design_matrix(train_design_block)
        print(f"    Z-scored design range: [{np.min(train_design_zscored.values):.4f}, {np.max(train_design_zscored.values):.4f}]")
        
        print(f"  Fitting GLM on training data...")
        print(f"    signal_scaling='psc' - applied manually")
        print(f"    standardize=True - applied manually to design matrix")
        print(f"    Run-specific nuisance regressors: {len(train_designs)} runs with independent drifts")
        
        fmri_glm = FirstLevelModel(
            t_r=tr,
            mask_img=mask_img,
            standardize=False,           # Already applied manually
            signal_scaling=False,         # Already applied manually
            hrf_model=hrf_model,
            minimize_memory=True,
            verbose=0
        )
        
        fmri_glm = fmri_glm.fit(train_img_scaled, design_matrices=train_design_zscored)
        
        # Extract task betas (only the first n_task_regressors columns)
        # Note: The task regressors are the same across runs, so we can use a simple contrast
        contrast_matrix = np.zeros((n_task_regressors, train_design_zscored.shape[1]))
        for i in range(n_task_regressors):
            contrast_matrix[i, i] = 1
        
        print(f"  Extracting task betas...")
        betas_img = fmri_glm.compute_contrast(
            contrast_matrix,
            output_type='effect_size'
        )
        betas = get_data(betas_img).astype(np.float32)
        print(f"    Betas shape: {betas.shape}")
        print(f"    Betas range: [{np.min(betas):.4f}, {np.max(betas):.4f}]")
        
        # --- TESTING ---
        print(f"\n  --- Testing Phase ---")
        test_imgs = [all_runs_imgs[i] for i in test_idx]
        test_designs = [designs[i] for i in test_idx]
        
        # Calculate total timepoints for testing
        test_total_tp = sum(run_timepoints[i] for i in test_idx)
        test_design_total_tp = sum(len(designs[i]) for i in test_idx)
        print(f"  Testing total timepoints: fMRI={test_total_tp}, Design={test_design_total_tp}")
        
        if test_total_tp != test_design_total_tp:
            raise ValueError(f"Testing timepoint mismatch! fMRI={test_total_tp}, Design={test_design_total_tp}")
        
        print(f"  Concatenating {len(test_imgs)} test runs...")
        test_imgs_concat = concat_imgs(test_imgs)
        test_data = get_data(test_imgs_concat)
        print(f"    Raw test data range: [{np.min(test_data):.2f}, {np.max(test_data):.2f}]")
        
        # Apply percent change scaling to test data
        print(f"  Applying percent change scaling to test data...")
        test_data_scaled = percent_change_scaling(test_data, mask_img)
        print(f"    Scaled test data range: [{np.min(test_data_scaled):.2f}, {np.max(test_data_scaled):.2f}]")
        
        # CRITICAL: Create block-diagonal design for test data too
        # This is needed for nuisance regression
        print(f"  Creating block-diagonal design matrix for test data...")
        test_design_block = create_block_diagonal_design(test_designs, n_task_regressors)
        print(f"    Test block design shape: {test_design_block.shape}")
        
        # Z-score the test design matrix (for nuisance regression)
        test_design_zscored = zscore_design_matrix(test_design_block)
        
        # Regress out nuisance regressors
        print(f"  Regressing out nuisance regressors...")
        
        nuisance_regressors = test_design_zscored.values[:, n_task_regressors:].astype(np.float32)
        print(f"    Nuisance regressors shape: {nuisance_regressors.shape}")
        print(f"    Number of nuisance regressors: {nuisance_regressors.shape[1]} (should be {len(test_designs)} * (original nuisance count))")
        
        confounds_df = pd.DataFrame(
            nuisance_regressors,
            columns=[f"nuisance_{i}" for i in range(nuisance_regressors.shape[1])]
        )
        
        # Create a Nifti image from scaled test data
        test_img_scaled = new_img_like(test_imgs_concat, test_data_scaled)
        
        # Use clean_img to remove nuisance effects
        test_imgs_cleaned = clean_img(
            test_img_scaled,
            confounds=confounds_df,
            detrend=False,                # Drift terms already in confounds
            standardize=False,             # Keep percent change scale
            t_r=tr
        )
        
        test_data_cleaned = get_data(test_imgs_cleaned).astype(np.float32)
        print(f"    Cleaned test data range: [{np.min(test_data_cleaned):.2f}, {np.max(test_data_cleaned):.2f}]")
        print(f"    Cleaned test data mean: {np.mean(test_data_cleaned):.2f}")
        print(f"    Cleaned test data std: {np.std(test_data_cleaned):.2f}")
        
        # Get task regressors for prediction (from z-scored design matrix)
        test_task_regressors = test_design_zscored.values[:, :n_task_regressors].astype(np.float32)
        print(f"    Test task regressors shape: {test_task_regressors.shape}")
        print(f"    Task regressors range: [{np.min(test_task_regressors):.4f}, {np.max(test_task_regressors):.4f}]")
        
        # Reshape for prediction
        betas_reshaped = betas.reshape(-1, n_task_regressors)
        mask_flat = mask_data.ravel()
        betas_masked = betas_reshaped[mask_flat, :]
        print(f"    Betas masked shape: {betas_masked.shape}")
        print(f"    Betas masked range: [{np.min(betas_masked):.4f}, {np.max(betas_masked):.4f}]")
        
        test_data_reshaped = test_data_cleaned.reshape(-1, test_data_cleaned.shape[-1])
        test_data_masked = test_data_reshaped[mask_flat, :]
        print(f"    Test data masked shape: {test_data_masked.shape}")
        print(f"    Test data masked range: [{np.min(test_data_masked):.4f}, {np.max(test_data_masked):.4f}]")
        print(f"    Test data masked mean: {np.mean(test_data_masked):.4f}")
        print(f"    Test data masked std: {np.std(test_data_masked):.4f}")
        
        # Predict
        print(f"  Predicting test time series...")
        predicted = np.zeros_like(test_data_masked)
        
        for task_idx in range(n_task_regressors):
            task_beta = betas_masked[:, task_idx:task_idx+1]
            task_regressor = test_task_regressors[:, task_idx:task_idx+1].T
            task_contribution = task_beta @ task_regressor
            predicted += task_contribution
            print(f"    Task {task_idx+1} contribution range: [{np.min(task_contribution):.4f}, {np.max(task_contribution):.4f}]")
        
        print(f"    Predicted range: [{np.min(predicted):.4f}, {np.max(predicted):.4f}]")
        print(f"    Predicted mean: {np.mean(predicted):.4f}")
        print(f"    Predicted std: {np.std(predicted):.4f}")
        
        # Compare scales
        scale_ratio = np.std(predicted) / (np.std(test_data_masked) + 1e-8)
        print(f"    Scale ratio (predicted/actual): {scale_ratio:.4f} (should be ~1)")
        
        # Compute R² with multiple methods
        print(f"  Computing R² with multiple methods...")
        
        test_data_4d = np.zeros(mask_shape + (test_data_masked.shape[1],), dtype=np.float32)
        predicted_4d = np.zeros(mask_shape + (predicted.shape[1],), dtype=np.float32)
        
        test_data_4d[mask_data, :] = test_data_masked
        predicted_4d[mask_data, :] = predicted
        
        r2_corr, r2_1minus, r2_sklearn = compute_r2_multiple_methods(
            test_data_4d, predicted_4d, mask_data
        )
        
        # Accumulate results for each method
        r2_corr_sum += r2_corr
        r2_1minus_sum += r2_1minus
        r2_sklearn_sum += r2_sklearn
        r2_corr_sq_sum += r2_corr ** 2

        # Clean up
        del fmri_glm, betas, betas_reshaped, betas_masked, predicted, test_data_cleaned
        del test_imgs_cleaned, train_imgs_concat, test_imgs_concat
        gc.collect()

    n_splits = len(splits)
    
    # Final results for each method
    mean_r2_corr = r2_corr_sum / n_splits
    mean_r2_1minus = r2_1minus_sum / n_splits
    mean_r2_sklearn = r2_sklearn_sum / n_splits
    var_r2_corr = (r2_corr_sq_sum / n_splits) - (mean_r2_corr ** 2)

    print(f"\n{'='*60}")
    print("FINAL RESULTS - ACROSS ALL SPLITS")
    print(f"{'='*60}")
    
    print(f"\nMethod 1 (Squared Correlation):")
    print(f"  Mean R²: {np.mean(mean_r2_corr[mask_data]):.6f}")
    print(f"  Std R²:  {np.std(mean_r2_corr[mask_data]):.6f}")
    print(f"  Range:   [{np.min(mean_r2_corr[mask_data]):.6f}, {np.max(mean_r2_corr[mask_data]):.6f}]")
    
    print(f"\nMethod 2 (1 - SS_res/SS_tot):")
    print(f"  Mean R²: {np.mean(mean_r2_1minus[mask_data]):.6f}")
    print(f"  Std R²:  {np.std(mean_r2_1minus[mask_data]):.6f}")
    print(f"  Range:   [{np.min(mean_r2_1minus[mask_data]):.6f}, {np.max(mean_r2_1minus[mask_data]):.6f}]")
    
    print(f"\nMethod 3 (sklearn r2_score):")
    print(f"  Mean R²: {np.mean(mean_r2_sklearn[mask_data]):.6f}")
    print(f"  Std R²:  {np.std(mean_r2_sklearn[mask_data]):.6f}")
    print(f"  Range:   [{np.min(mean_r2_sklearn[mask_data]):.6f}, {np.max(mean_r2_sklearn[mask_data]):.6f}]")
    
    # Cross-method correlations for final results
    corr_1_2 = np.corrcoef(mean_r2_corr[mask_data].ravel(), mean_r2_1minus[mask_data].ravel())[0,1]
    corr_1_3 = np.corrcoef(mean_r2_corr[mask_data].ravel(), mean_r2_sklearn[mask_data].ravel())[0,1]
    print(f"\nCross-method correlations (final maps):")
    print(f"  Method 1 vs 2: {corr_1_2:.6f}")
    print(f"  Method 1 vs 3: {corr_1_3:.6f}")

    return mean_r2_corr.astype(np.float32), var_r2_corr.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Cross-validated GLM with block-diagonal design matrices")
    
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
    print("CROSS-VALIDATED GLM WITH BLOCK-DIAGONAL DESIGN MATRICES")
    print(f"{'='*60}")
    print(f"Input files: {len(args.nii_files)} runs")
    print(f"Parameters:")
    print(f"  Task regressors: {args.n_task_regressors}")
    print(f"  TR: {args.tr}s")
    print(f"  HRF model: {args.hrf_model}")
    print(f"  Number of splits: {args.n_splits}")
    print(f"  Test size: {args.test_size*100:.0f}%")
    print(f"  Manual percent change scaling applied to all data")
    print(f"  Manual z-scoring applied to design matrices")
    print(f"  Block-diagonal design with run-specific nuisance regressors")
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
        print
