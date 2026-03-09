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
from tqdm import tqdm
import gc
import warnings
warnings.filterwarnings('ignore')


# ============================================================
# Utilities
# ============================================================

def percent_change_scaling_per_run(img, mask_img, debug=False):
    """
    Apply percent change scaling PER RUN: (x - mean)/mean * 100
    This is done independently for each run before concatenation
    """
    # Get data as numpy array
    data = get_data(img).astype(np.float64)
    mask_data = mask_img.get_fdata().astype(bool)
    
    # Get original shape
    original_shape = data.shape
    n_timepoints = original_shape[-1]
    
    if debug:
        print(f"\n    DEBUG - Percent Change Scaling (per run):")
        print(f"      Original data shape: {original_shape}")
        print(f"      Original data range: [{np.min(data):.2f}, {np.max(data):.2f}]")
        print(f"      Original data mean: {np.mean(data):.2f}")
    
    # Reshape to voxels × time
    n_voxels = np.prod(original_shape[:-1])
    data_reshaped = data.reshape(n_voxels, n_timepoints)
    mask_flat = mask_data.ravel()
    
    # Get masked data
    masked_data = data_reshaped[mask_flat, :]
    
    if debug and masked_data.size > 0:
        print(f"      Masked data shape: {masked_data.shape}")
        print(f"      Masked data range: [{np.min(masked_data):.2f}, {np.max(masked_data):.2f}]")
        print(f"      Masked data mean: {np.mean(masked_data):.2f}")
    
    if masked_data.size > 0:
        # Compute mean for each voxel WITHIN THIS RUN ONLY
        voxel_means = np.mean(masked_data, axis=1, keepdims=True)
        
        if debug:
            print(f"      Voxel means range: [{np.min(voxel_means):.2f}, {np.max(voxel_means):.2f}]")
            print(f"      Voxel means mean: {np.mean(voxel_means):.2f}")
        
        # Avoid division by zero
        voxel_means = np.where(np.abs(voxel_means) < 1e-6, 1, voxel_means)
        
        # Apply percent change scaling to masked voxels
        masked_pc = (masked_data - voxel_means) / voxel_means * 100
        
        if debug:
            print(f"      Percent change range: [{np.min(masked_pc):.2f}, {np.max(masked_pc):.2f}]")
            print(f"      Percent change mean: {np.mean(masked_pc):.2f}")
            print(f"      Percent change std: {np.std(masked_pc):.2f}")
        
        # Put back into full array
        data_pc_reshaped = np.zeros_like(data_reshaped)
        data_pc_reshaped[mask_flat, :] = masked_pc
    else:
        data_pc_reshaped = data_reshaped
    
    # Reshape back to original 4D
    data_pc = data_pc_reshaped.reshape(original_shape).astype(np.float32)
    
    return new_img_like(img, data_pc)


def load_and_preprocess_runs(nii_files, mask_img, debug=False):
    """
    Load all runs and apply percent change scaling PER RUN
    This is critical - each run is scaled independently before any concatenation
    """
    print("\n" + "="*60)
    print("STEP 1: Loading runs and applying per-run percent change scaling")
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
        
        # CRITICAL: Apply percent change scaling PER RUN, before any concatenation
        print(f"    Applying percent change scaling (per run)...")
        img_pc = percent_change_scaling_per_run(img, mask_img, debug=debug)
        all_runs_imgs.append(img_pc)
        
        print(f"    Shape: {img_pc.shape}, Timepoints: {n_timepoints}")
        
        # Quick check of scaled data
        if debug:
            data_sample = get_data(img_pc)[mask_img.get_fdata().astype(bool)]
            if len(data_sample) > 0:
                print(f"    Scaled data stats (within mask):")
                print(f"      Range: [{np.min(data_sample):.2f}, {np.max(data_sample):.2f}]")
                print(f"      Mean: {np.mean(data_sample):.2f}")
                print(f"      Std: {np.std(data_sample):.2f}")
    
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
            print(f"    Task regressor range: [{np.min(df.values[:, :n_task_regressors]):.4f}, {np.max(df.values[:, :n_task_regressors]):.4f}]")

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
    
    # Initialize arrays for each method
    r2_corr = np.zeros(n_voxels)
    r2_1minus = np.zeros(n_voxels)
    r2_sklearn = np.zeros(n_voxels)
    
    # Additional diagnostics
    neg_count_1minus = 0
    neg_count_sklearn = 0
    
    for v in range(n_voxels):
        y_t = y_true_valid[v, :]
        y_p = y_pred_valid[v, :]
        
        # Method 1: Squared Pearson correlation
        corr = np.corrcoef(y_t, y_p)[0, 1]
        r2_corr[v] = corr ** 2 if not np.isnan(corr) else 0
        
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
    print(f"      Correlation between Method 1 and 2: {np.corrcoef(r2_corr, r2_1minus)[0,1]:.6f}")
    print(f"      Correlation between Method 1 and 3: {np.corrcoef(r2_corr, r2_sklearn)[0,1]:.6f}")
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
    debug=False,
):
    """
    Run cross-validated GLM with per-run percent change scaling
    """
    
    # Step 1: Load all runs and apply per-run percent change scaling
    all_runs_imgs, reference_img, run_timepoints = load_and_preprocess_runs(nii_files, mask_img, debug=debug)
    
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
    
    # Store results
    r2_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sq_sum = np.zeros(mask_shape, dtype=np.float32)

    for split_id, (train_idx, test_idx) in enumerate(splits):
        
        print(f"\n{'='*60}")
        print(f"STEP 3: Processing split {split_id + 1}/{len(splits)}")
        print(f"{'='*60}")
        print(f"  Training runs: {train_idx}")
        print(f"  Testing runs: {test_idx}")
        
        # --- TRAINING ---
        print(f"\n  --- Training Phase ---")
        # Get training images (ALREADY in percent change from per-run scaling)
        train_imgs = [all_runs_imgs[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        
        # Calculate total timepoints for training
        train_total_tp = sum(run_timepoints[i] for i in train_idx)
        train_design_total_tp = sum(len(designs[i]) for i in train_idx)
        print(f"  Training total timepoints: fMRI={train_total_tp}, Design={train_design_total_tp}")
        
        if train_total_tp != train_design_total_tp:
            raise ValueError(f"Training timepoint mismatch! fMRI={train_total_tp}, Design={train_design_total_tp}")
        
        # Concatenate training images (already scaled per run)
        print(f"  Concatenating {len(train_imgs)} training runs...")
        train_imgs_concat = concat_imgs(train_imgs)
        train_data = get_data(train_imgs_concat)
        print(f"    Concatenated training data shape: {train_data.shape}")
        print(f"    Concatenated training data range: [{np.min(train_data):.2f}, {np.max(train_data):.2f}]")
        print(f"    Concatenated training data mean: {np.mean(train_data):.2f}")
        
        # Create proper block-diagonal design matrix with run-specific nuisance regressors
        print(f"  Creating block-diagonal design matrix...")
        
        # For training, we need to stack the design matrices with run-specific nuisance regressors
        train_design_blocks = []
        for run_idx, design in enumerate(train_designs):
            # For each run, keep task regressors as is, but make nuisance regressors run-specific
            # by adding a suffix to their names
            task_part = design.iloc[:, :n_task_regressors]
            nuisance_part = design.iloc[:, n_task_regressors:]
            
            # Rename nuisance columns to be run-specific
            rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
            nuisance_part = nuisance_part.rename(columns=rename_dict)
            
            # Combine task and renamed nuisance
            run_design = pd.concat([task_part, nuisance_part], axis=1)
            train_design_blocks.append(run_design)
        
        # Concatenate all runs
        train_design_matrix = pd.concat(train_design_blocks, axis=0, ignore_index=True)
        train_design_matrix = train_design_matrix.fillna(0)
        
        print(f"    Training design matrix shape: {train_design_matrix.shape}")
        print(f"    Number of regressors: {train_design_matrix.shape[1]}")
        print(f"    Training design matrix range: [{np.min(train_design_matrix.values):.4f}, {np.max(train_design_matrix.values):.4f}]")
        
        # Fit GLM with signal_scaling=False since we already applied percent change
        print(f"  Fitting GLM on training data...")
        print(f"    signal_scaling=False - data already in percent change")
        print(f"    standardize=True - z-scores design matrix")
        
        fmri_glm = FirstLevelModel(
            t_r=tr,
            mask_img=mask_img,
            standardize=True,
            signal_scaling=False,  # CRITICAL: Data already in percent change
            hrf_model=hrf_model,
            minimize_memory=True,
            verbose=0
        )
        
        fmri_glm = fmri_glm.fit(train_imgs_concat, design_matrices=train_design_matrix)
        
        # Extract task betas (first n_task_regressors)
        contrast_matrix = np.zeros((n_task_regressors, train_design_matrix.shape[1]))
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
        print(f"    Betas mean: {np.mean(betas):.4f}")
        
        # Store training task regressor statistics for later z-scoring of test data
        train_task_means = []
        train_task_stds = []
        for i in range(n_task_regressors):
            train_task_col = train_design_matrix.values[:, i]
            train_task_means.append(np.mean(train_task_col))
            train_task_stds.append(np.std(train_task_col))
        
        print(f"    Training task regressor stats:")
        for i in range(n_task_regressors):
            print(f"      Task {i}: mean={train_task_means[i]:.4f}, std={train_task_stds[i]:.4f}")
        
        # --- TESTING ---
        print(f"\n  --- Testing Phase ---")
        # Get test images (ALREADY in percent change from per-run scaling)
        test_imgs = [all_runs_imgs[i] for i in test_idx]
        test_designs = [designs[i] for i in test_idx]
        
        # Calculate total timepoints for testing
        test_total_tp = sum(run_timepoints[i] for i in test_idx)
        test_design_total_tp = sum(len(designs[i]) for i in test_idx)
        print(f"  Testing total timepoints: fMRI={test_total_tp}, Design={test_design_total_tp}")
        
        if test_total_tp != test_design_total_tp:
            raise ValueError(f"Testing timepoint mismatch! fMRI={test_total_tp}, Design={test_design_total_tp}")
        
        # Concatenate test images (already scaled per run)
        print(f"  Concatenating {len(test_imgs)} test runs...")
        test_imgs_concat = concat_imgs(test_imgs)
        test_data = get_data(test_imgs_concat)
        print(f"    Concatenated test data shape: {test_data.shape}")
        print(f"    Concatenated test data range: [{np.min(test_data):.2f}, {np.max(test_data):.2f}]")
        print(f"    Concatenated test data mean: {np.mean(test_data):.2f}")
        
        # Create test design matrix with run-specific nuisance regressors
        test_design_blocks = []
        for run_idx, design in enumerate(test_designs):
            task_part = design.iloc[:, :n_task_regressors]
            nuisance_part = design.iloc[:, n_task_regressors:]
            
            # Rename nuisance columns to be run-specific
            rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
            nuisance_part = nuisance_part.rename(columns=rename_dict)
            
            # Combine task and renamed nuisance
            run_design = pd.concat([task_part, nuisance_part], axis=1)
            test_design_blocks.append(run_design)
        
        test_design_matrix = pd.concat(test_design_blocks, axis=0, ignore_index=True)
        test_design_matrix = test_design_matrix.fillna(0)
        print(f"    Test design matrix shape: {test_design_matrix.shape}")
        print(f"    Test design matrix range: [{np.min(test_design_matrix.values):.4f}, {np.max(test_design_matrix.values):.4f}]")
        
        # Extract nuisance regressors for cleaning
        nuisance_regressors = test_design_matrix.values[:, n_task_regressors:].astype(np.float32)
        print(f"    Nuisance regressors shape: {nuisance_regressors.shape}")
        print(f"    Nuisance regressors range: [{np.min(nuisance_regressors):.4f}, {np.max(nuisance_regressors):.4f}]")
        
        # Create confounds DataFrame
        confounds_df = pd.DataFrame(
            nuisance_regressors,
            columns=[f"nuisance_{i}" for i in range(nuisance_regressors.shape[1])]
        )
        
        # Clean test data by regressing out nuisance
        print(f"  Regressing out nuisance regressors...")
        test_imgs_cleaned = clean_img(
            test_imgs_concat,  # Use already scaled data
            confounds=confounds_df,
            detrend=False,
            standardize=False,  # Keep percent change scale
            t_r=tr
        )
        
        test_data_cleaned = get_data(test_imgs_cleaned).astype(np.float32)
        print(f"    Cleaned test data range: [{np.min(test_data_cleaned):.2f}, {np.max(test_data_cleaned):.2f}]")
        print(f"    Cleaned test data mean: {np.mean(test_data_cleaned):.2f}")
        print(f"    Cleaned test data std: {np.std(test_data_cleaned):.2f}")
        
        # Check percentiles to identify outliers
        test_data_masked_temp = test_data_cleaned[mask_data]
        if len(test_data_masked_temp) > 0:
            percentiles = np.percentile(test_data_masked_temp, [1, 5, 50, 95, 99])
            print(f"    Cleaned data percentiles: 1%={percentiles[0]:.2f}, 5%={percentiles[1]:.2f}, "
                  f"50%={percentiles[2]:.2f}, 95%={percentiles[3]:.2f}, 99%={percentiles[4]:.2f}")
        
        # Get raw task regressors from test design
        test_task_regressors_raw = test_design_matrix.values[:, :n_task_regressors].astype(np.float32)
        
        # CRITICAL: Z-score test task regressors using training statistics
        # This is necessary because the GLM was fit with standardize=True
        test_task_regressors_zscored = np.zeros_like(test_task_regressors_raw)
        for i in range(n_task_regressors):
            if train_task_stds[i] > 0:
                test_task_regressors_zscored[:, i] = (test_task_regressors_raw[:, i] - train_task_means[i]) / train_task_stds[i]
            else:
                test_task_regressors_zscored[:, i] = 0
        
        print(f"    Raw task regressors range: [{np.min(test_task_regressors_raw):.4f}, {np.max(test_task_regressors_raw):.4f}]")
        print(f"    Z-scored task regressors range: [{np.min(test_task_regressors_zscored):.4f}, {np.max(test_task_regressors_zscored):.4f}]")
        
        # Reshape for prediction
        betas_reshaped = betas.reshape(-1, n_task_regressors)
        mask_flat = mask_data.ravel()
        betas_masked = betas_reshaped[mask_flat, :]
        
        test_data_reshaped = test_data_cleaned.reshape(-1, test_data_cleaned.shape[-1])
        test_data_masked = test_data_reshaped[mask_flat, :]
        
        # Predict using Z-SCORED task regressors
        print(f"  Predicting test time series with z-scored task regressors...")
        predicted = betas_masked @ test_task_regressors_zscored.T
        
        print(f"    Predicted range: [{np.min(predicted):.4f}, {np.max(predicted):.4f}]")
        print(f"    Predicted mean: {np.mean(predicted):.4f}")
        print(f"    Predicted std: {np.std(predicted):.4f}")
        print(f"    Cleaned data range: [{np.min(test_data_masked):.4f}, {np.max(test_data_masked):.4f}]")
        print(f"    Cleaned data mean: {np.mean(test_data_masked):.4f}")
        print(f"    Cleaned data std: {np.std(test_data_masked):.4f}")
        
        # Scale comparison
        scale_ratio = np.std(predicted) / (np.std(test_data_masked) + 1e-8)
        print(f"    Scale ratio (predicted/actual): {scale_ratio:.4f}")
        
        # Compute R²
        print(f"  Computing R²...")
        
        test_data_4d = np.zeros(mask_shape + (test_data_masked.shape[1],), dtype=np.float32)
        predicted_4d = np.zeros(mask_shape + (predicted.shape[1],), dtype=np.float32)
        
        test_data_4d[mask_data, :] = test_data_masked
        predicted_4d[mask_data, :] = predicted
        
        r2_map, r2_1minus, r2_sklearn = compute_r2_multiple_methods(test_data_4d, predicted_4d, mask_data)
        
        # Accumulate results (using squared correlation as primary)
        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2

        # Clean up
        del fmri_glm, betas, betas_reshaped, betas_masked, predicted, test_data_cleaned
        del test_imgs_cleaned, train_imgs_concat, test_imgs_concat
        gc.collect()

    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)

    print(f"\n{'='*60}")
    print("FINAL RESULTS")
    print(f"{'='*60}")
    print(f"Mean R² across all splits (within mask): {np.mean(mean_r2[mask_data]):.6f}")
    print(f"Std R² across all splits: {np.std(mean_r2[mask_data]):.6f}")
    print(f"R² range: [{np.min(mean_r2[mask_data]):.6f}, {np.max(mean_r2[mask_data]):.6f}]")

    return mean_r2.astype(np.float32), var_r2.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Cross-validated GLM with per-run percent change scaling")
    
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
    parser.add_argument("--debug", action="store_true",
                        help="Enable debug output")

    args = parser.parse_args()

    if len(args.nii_files) != len(args.events_files):
        raise ValueError(f"Number of NIfTI files ({len(args.nii_files)}) does not match number of event files ({len(args.events_files)})")

    print(f"\n{'='*60}")
    print("CROSS-VALIDATED GLM WITH PER-RUN PERCENT CHANGE SCALING")
    print(f"{'='*60}")
    print(f"Input files: {len(args.nii_files)} runs")
    print(f"Parameters:")
    print(f"  Task regressors: {args.n_task_regressors}")
    print(f"  TR: {args.tr}s")
    print(f"  HRF model: {args.hrf_model}")
    print(f"  Number of splits: {args.n_splits}")
    print(f"  Test size: {args.test_size*100:.0f}%")
    print(f"  Debug mode: {args.debug}")
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
            permute=False,
            debug=args.debug,
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
                permute=True,
                debug=args.debug,
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

