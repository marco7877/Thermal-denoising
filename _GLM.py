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
        print(f"      Original data mean (all voxels+time): {np.mean(data):.2f}")
    
    # Reshape to voxels × time
    n_voxels = np.prod(original_shape[:-1])
    data_reshaped = data.reshape(n_voxels, n_timepoints)
    mask_flat = mask_data.ravel()
    
    # Get masked data
    masked_data = data_reshaped[mask_flat, :].copy()  # Make a copy to avoid reference issues
    
    if debug and masked_data.size > 0:
        print(f"      Masked data shape: {masked_data.shape}")
        print(f"      Masked data range (all voxels+time): [{np.min(masked_data):.2f}, {np.max(masked_data):.2f}]")
        print(f"      Masked data mean (all voxels+time): {np.mean(masked_data):.2f}")
    
    if masked_data.size > 0:
        # Compute mean for each voxel WITHIN THIS RUN ONLY
        voxel_means = np.mean(masked_data, axis=1, keepdims=True)
        
        if debug:
            print(f"      Voxel means range (per voxel baseline): [{np.min(voxel_means):.2f}, {np.max(voxel_means):.2f}]")
            print(f"      Average of voxel means: {np.mean(voxel_means):.2f}")
            print(f"      Mean per voxel shape: {voxel_means.shape}")
        
        # Avoid division by zero
        voxel_means_safe = np.where(np.abs(voxel_means) < 1e-6, 1, voxel_means)
        
        # Apply percent change scaling to masked voxels
        # This is the key line - we need to broadcast correctly
        masked_pc = (masked_data - voxel_means_safe) / voxel_means_safe * 100
        
        # Calculate per-voxel means after transformation
        per_voxel_pc_means = np.mean(masked_pc, axis=1)
        
        if debug:
            print(f"      Percent change range (all voxels+time): [{np.min(masked_pc):.2f}, {np.max(masked_pc):.2f}]")
            print(f"      Percent change mean (all voxels+time): {np.mean(masked_pc):.2f}")
            print(f"      Per-voxel percent change means - range: [{np.min(per_voxel_pc_means):.2f}, {np.max(per_voxel_pc_means):.2f}]")
            print(f"      Per-voxel percent change means - average: {np.mean(per_voxel_pc_means):.2f} (should be 0)")
            
            # Check if all per-voxel means are the same (which would indicate a bug)
            if np.std(per_voxel_pc_means) < 1e-6:
                print(f"      ERROR: All per-voxel means are identical! This indicates a broadcasting bug.")
                print(f"      First 5 per-voxel means: {per_voxel_pc_means[:5]}")
            
            # This is the key check - the average of per-voxel means should be 0
            if abs(np.mean(per_voxel_pc_means)) > 0.01:
                print(f"      ERROR: Per-voxel means average is {np.mean(per_voxel_pc_means):.2f}, should be 0!")
                print(f"      This indicates a bug in the percent change calculation")
                
                # Let's debug by looking at a single voxel
                test_voxel = 0
                print(f"      Debug single voxel {test_voxel}:")
                print(f"        Original data (first 5): {masked_data[test_voxel, :5]}")
                print(f"        Voxel mean: {voxel_means[test_voxel, 0]:.2f}")
                pc_test = (masked_data[test_voxel, :] - voxel_means[test_voxel, 0]) / voxel_means[test_voxel, 0] * 100
                print(f"        Calculated PC (first 5): {pc_test[:5]}")
                print(f"        Mean of PC for this voxel: {np.mean(pc_test):.4f}")
        
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
# Single split function
# ============================================================

def run_single_split(train_imgs, test_imgs, train_designs, test_designs,
                     mask_data, n_task_regressors, tr, hrf_model):
    """
    Run one CV split and return the R² map.
    """
    # Concatenate training images
    train_imgs_concat = concat_imgs(train_imgs)
    
    # Create block-diagonal design matrix for training
    train_design_blocks = []
    for run_idx, design in enumerate(train_designs):
        task_part = design.iloc[:, :n_task_regressors]
        nuisance_part = design.iloc[:, n_task_regressors:]
        rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
        nuisance_part = nuisance_part.rename(columns=rename_dict)
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        train_design_blocks.append(run_design)
    train_design_matrix = pd.concat(train_design_blocks, axis=0, ignore_index=True).fillna(0)
    
    # Fit GLM
    fmri_glm = FirstLevelModel(
        t_r=tr,
        mask_img=None,
        standardize=True,
        signal_scaling=False,
        hrf_model=hrf_model,
        minimize_memory=True
    )
    fmri_glm = fmri_glm.fit(train_imgs_concat, design_matrices=train_design_matrix)
    
    # Extract task betas
    contrast_matrix = np.zeros((n_task_regressors, train_design_matrix.shape[1]))
    for i in range(n_task_regressors):
        contrast_matrix[i, i] = 1
    betas_img = fmri_glm.compute_contrast(contrast_matrix, output_type='effect_size')
    betas = get_data(betas_img).astype(np.float32)
    
    # Training task stats for z-scoring test regressors
    train_task_means = [np.mean(train_design_matrix.values[:, i]) for i in range(n_task_regressors)]
    train_task_stds = [np.std(train_design_matrix.values[:, i]) for i in range(n_task_regressors)]
    
    # --- Testing (using real test designs) ---
    test_imgs_concat = concat_imgs(test_imgs)
    
    # Create test design matrix (nuisance renamed, but task part is real)
    test_design_blocks = []
    for run_idx, design in enumerate(test_designs):
        task_part = design.iloc[:, :n_task_regressors]
        nuisance_part = design.iloc[:, n_task_regressors:]
        rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
        nuisance_part = nuisance_part.rename(columns=rename_dict)
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        test_design_blocks.append(run_design)
    test_design_matrix = pd.concat(test_design_blocks, axis=0, ignore_index=True).fillna(0)
    
    # Regress out nuisance from test data
    nuisance_regressors = test_design_matrix.values[:, n_task_regressors:].astype(np.float32)
    confounds_df = pd.DataFrame(nuisance_regressors, columns=[f"nuisance_{i}" for i in range(nuisance_regressors.shape[1])])
    test_imgs_cleaned = clean_img(test_imgs_concat, confounds=confounds_df, detrend=False, standardize=False, t_r=tr)
    test_data_cleaned = get_data(test_imgs_cleaned).astype(np.float32)
    
    # Z-score test task regressors using training stats
    test_task_regressors_raw = test_design_matrix.values[:, :n_task_regressors].astype(np.float32)
    test_task_regressors_zscored = np.zeros_like(test_task_regressors_raw)
    for i in range(n_task_regressors):
        if train_task_stds[i] > 0:
            test_task_regressors_zscored[:, i] = (test_task_regressors_raw[:, i] - train_task_means[i]) / train_task_stds[i]
    
    # Predict
    betas_reshaped = betas.reshape(-1, n_task_regressors)
    mask_flat = mask_data.ravel()
    betas_masked = betas_reshaped[mask_flat, :]
    test_data_reshaped = test_data_cleaned.reshape(-1, test_data_cleaned.shape[-1])
    test_data_masked = test_data_reshaped[mask_flat, :]
    predicted = betas_masked @ test_task_regressors_zscored.T
    
    # Compute R² map
    r2_map = np.zeros(mask_data.shape, dtype=np.float32)
    for v in range(test_data_masked.shape[0]):
        y_t = test_data_masked[v, :]
        y_p = predicted[v, :]
        corr = np.corrcoef(y_t, y_p)[0, 1]
        r2_map[mask_data][v] = corr ** 2 if not np.isnan(corr) else 0
    
    return r2_map


# ============================================================
# Optimized permutation function (permute training designs only)
# ============================================================

def run_permutations(
    nii_files,
    events_files,
    mask_img,
    n_task_regressors,
    tr,
    hrf_model,
    splits,
    n_permutations,
    random_state,
    debug=False,
):
    """
    Run permutations by shuffling only the training design matrices.
    Data is preprocessed once and reused.
    """
    # Preprocess all data once
    all_runs_imgs, reference_img, run_timepoints = load_and_preprocess_runs(nii_files, mask_img, debug=debug)
    designs = load_design_matrices(events_files, run_timepoints, n_task_regressors)
    mask_data = mask_img.get_fdata().astype(bool)
    
    max_distribution = np.zeros(n_permutations, dtype=np.float32)
    
    for p in tqdm(range(n_permutations), desc="Permutations"):
        rng = np.random.default_rng(random_state + p + 1)
        r2_perm_sum = np.zeros(mask_data.shape, dtype=np.float32)
         for split_id, (train_idx, test_idx) in enumerate(splits):
            # Test data (unchanged)
            print(f"Split {split_id+1}/{len(splits)}")
            test_imgs = [all_runs_imgs[i] for i in test_idx]
            test_designs = [designs[i] for i in test_idx]
            
            # Training data (images unchanged, designs permuted)
            train_imgs = [all_runs_imgs[i] for i in train_idx]
            train_designs_orig = [designs[i] for i in train_idx]
            
            # Permute rows of each training design matrix
            train_designs_perm = []
            for design in train_designs_orig:
                design_perm = design.sample(frac=1, random_state=rng).reset_index(drop=True)
                train_designs_perm.append(design_perm)
            
            # Run the split
            r2_map = run_single_split(
                train_imgs, test_imgs,
                train_designs_perm, test_designs,
                mask_data, n_task_regressors, tr, hrf_model
            )
            r2_perm_sum += r2_map
        
        # Average over splits (if multiple splits per permutation)
        r2_perm_mean = r2_perm_sum / len(splits)
        max_distribution[p] = np.max(r2_perm_mean[mask_data])
        
        # Clean up
        del r2_perm_sum
        gc.collect()
    
    return max_distribution, r2_perm_mean.astype(np.float32)


# ============================================================
# Real data CV (unchanged)
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
    debug=False,
):
    """
    Run cross-validation for real data (no permutation).
    """
    all_runs_imgs, reference_img, run_timepoints = load_and_preprocess_runs(nii_files, mask_img, debug=debug)
    designs = load_design_matrices(events_files, run_timepoints, n_task_regressors)
    mask_data = mask_img.get_fdata().astype(bool)
    mask_shape = mask_data.shape
    
    r2_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sq_sum = np.zeros(mask_shape, dtype=np.float32)
    
    for split_id, (train_idx, test_idx) in enumerate(splits):
        print(f"Split {split_id + 1}/{len(splits)}")
        
        train_imgs = [all_runs_imgs[i] for i in train_idx]
        train_designs = [designs[i] for i in train_idx]
        test_imgs = [all_runs_imgs[i] for i in test_idx]
        test_designs = [designs[i] for i in test_idx]
        
        r2_map = run_single_split(
            train_imgs, test_imgs,
            train_designs, test_designs,
            mask_data, n_task_regressors, tr, hrf_model
        )
        
        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2
        gc.collect()
    
    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)
    
    print("FINAL RESULTS")
    print(f"Mean R²: {np.mean(mean_r2[mask_data]):.6f}")
    print(f"Std R²: {np.std(mean_r2[mask_data]):.6f}")
    print(f"Range: [{np.min(mean_r2[mask_data]):.6f}, {np.max(mean_r2[mask_data]):.6f}]")
    
    return mean_r2.astype(np.float32), var_r2.astype(np.float32)


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Cross-validated GLM")
    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_task_regressors", type=int, required=True)
    parser.add_argument("--tr", type=float, default=2.0)
    parser.add_argument("--hrf_model", type=str, default="spm")
    parser.add_argument("--n_splits", type=int, default=50)
    parser.add_argument("--test_size", type=float, default=0.3)
    parser.add_argument("--n_permutations", type=int, default=0)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--only_permutations", action="store_true")
    parser.add_argument("--debug", action="store_true")
    
    args = parser.parse_args()
    
    if len(args.nii_files) != len(args.events_files):
        raise ValueError("Mismatch between NIfTI and event files")
    
    print(f"\n{'='*60}")
    print("CROSS-VALIDATED GLM")
    print(f"{'='*60}")
    # ... print parameters
    
    mask_img = nib.load(args.mask)
    
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))
    
    if not args.only_permutations:
        print("\nRUNNING REAL DATA CV")
        mean_r2, var_r2 = run_cv(
            args.nii_files, args.events_files, mask_img,
            args.n_task_regressors, args.tr, args.hrf_model,
            splits, args.random_state, debug=args.debug
        )
        nib.Nifti1Image(mean_r2, mask_img.affine).to_filename(f"{args.output_prefix}_real_mean_r2.nii.gz")
        nib.Nifti1Image(var_r2, mask_img.affine).to_filename(f"{args.output_prefix}_real_var_r2.nii.gz")
        del mean_r2, var_r2
        gc.collect()
    
    if args.n_permutations > 0:
        print(f"\nRUNNING {args.n_permutations} PERMUTATIONS")
        max_dist, mean_r2 = run_permutations(
            args.nii_files, args.events_files, mask_img,
            args.n_task_regressors, args.tr, args.hrf_model,
            splits, args.n_permutations, args.random_state,
            debug=args.debug
        )
        np.save(f"{args.output_prefix}_perm_max_distribution.npy", max_dist)
        nib.Nifti1Image(mean_r2, mask_img.affine).to_filename(f"{args.output_prefix}_perm_mean_r2.nii.gz")
        del mean_r2, max_dist
        gc.collect()

    print("\nDONE!\n")


if __name__ == "__main__":
    main()
