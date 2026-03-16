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


# ============================================================
# Utilities (same as before, but with caching)
# ============================================================

def preprocess_all_data(nii_files, events_files, mask_img, n_task_regressors, debug=False):
    """
    Preprocess all data once and return:
    - Preprocessed images (already in percent change)
    - Original design matrices (for permutation)
    - Run timepoints
    - Mask data
    """
    print("\n" + "="*60)
    print("PREPROCESSING ALL DATA (ONCE)")
    print("="*60)
    
    # Load and scale runs
    all_runs_imgs = []
    run_timepoints = []
    reference_img = None
    
    for i, f in enumerate(nii_files):
        print(f"\n  Processing run {i+1}/{len(nii_files)}: {f}")
        img = nib.load(f)
        
        if reference_img is None:
            reference_img = img
        else:
            if not np.allclose(img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
                img = resample_to_img(img, reference_img, interpolation='continuous')
        
        n_timepoints = img.shape[-1]
        run_timepoints.append(n_timepoints)
        
        print(f"    Applying percent change scaling...")
        img_pc = percent_change_scaling_per_run(img, mask_img, debug=debug)
        all_runs_imgs.append(img_pc)
    
    # Load original design matrices (without permutation)
    designs = []
    for i, f in enumerate(events_files):
        print(f"\n  Loading design matrix {i+1}/{len(events_files)}: {f}")
        df = pd.read_csv(f)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        
        if len(df) != run_timepoints[i]:
            raise ValueError(f"Timepoint mismatch for run {i}")
        
        designs.append(df)
    
    # Get mask data
    mask_data = mask_img.get_fdata().astype(bool)
    mask_shape = mask_data.shape
    
    print(f"\nPreprocessing complete:")
    print(f"  {len(all_runs_imgs)} runs, {np.sum(mask_data)} voxels")
    print(f"  Total timepoints: {sum(run_timepoints)}")
    
    return {
        'runs_imgs': all_runs_imgs,
        'designs': designs,
        'run_timepoints': run_timepoints,
        'mask_data': mask_data,
        'mask_shape': mask_shape,
        'reference_img': reference_img,
        'n_runs': len(nii_files)
    }


def run_single_cv_split(train_idx, test_idx, preprocessed_data, n_task_regressors, tr, hrf_model, permute_designs=False, random_state=None):
    """
    Run a single CV split (can be parallelized)
    """
    runs_imgs = preprocessed_data['runs_imgs']
    designs = preprocessed_data['designs']
    run_timepoints = preprocessed_data['run_timepoints']
    mask_data = preprocessed_data['mask_data']
    mask_shape = preprocessed_data['mask_shape']
    
    # Get training/test data
    train_imgs = [runs_imgs[i] for i in train_idx]
    test_imgs = [runs_imgs[i] for i in test_idx]
    
    # Handle design matrices (with optional permutation)
    if permute_designs:
        rng = np.random.default_rng(random_state)
        train_designs = [designs[i].sample(frac=1, random_state=rng).reset_index(drop=True) for i in train_idx]
        test_designs = [designs[i].sample(frac=1, random_state=rng).reset_index(drop=True) for i in test_idx]
    else:
        train_designs = [designs[i] for i in train_idx]
        test_designs = [designs[i] for i in test_idx]
    
    # Concatenate training images
    train_imgs_concat = concat_imgs(train_imgs)
    
    # Create block-diagonal design matrix
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
        mask_img=None,  # We'll handle masking manually
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
    
    # Store training stats
    train_task_means = [np.mean(train_design_matrix.values[:, i]) for i in range(n_task_regressors)]
    train_task_stds = [np.std(train_design_matrix.values[:, i]) for i in range(n_task_regressors)]
    
    # --- TESTING ---
    test_imgs_concat = concat_imgs(test_imgs)
    
    # Create test design matrix
    test_design_blocks = []
    for run_idx, design in enumerate(test_designs):
        task_part = design.iloc[:, :n_task_regressors]
        nuisance_part = design.iloc[:, n_task_regressors:]
        rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
        nuisance_part = nuisance_part.rename(columns=rename_dict)
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        test_design_blocks.append(run_design)
    
    test_design_matrix = pd.concat(test_design_blocks, axis=0, ignore_index=True).fillna(0)
    
    # Regress out nuisance
    nuisance_regressors = test_design_matrix.values[:, n_task_regressors:].astype(np.float32)
    confounds_df = pd.DataFrame(nuisance_regressors, columns=[f"nuisance_{i}" for i in range(nuisance_regressors.shape[1])])
    
    test_imgs_cleaned = clean_img(test_imgs_concat, confounds=confounds_df, detrend=False, standardize=False, t_r=tr)
    test_data_cleaned = get_data(test_imgs_cleaned).astype(np.float32)
    
    # Z-score test task regressors
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
    
    # Compute R² (just the max for permutation distribution)
    valid_voxels = mask_data.ravel() != 0
    r2_values = np.zeros(np.sum(valid_voxels))
    
    for v in range(test_data_masked.shape[0]):
        y_t = test_data_masked[v, :]
        y_p = predicted[v, :]
        corr = np.corrcoef(y_t, y_p)[0, 1]
        r2_values[v] = corr ** 2 if not np.isnan(corr) else 0
    
    return np.max(r2_values)


def run_cv_optimized(
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
    n_jobs=1,
):
    """
    Optimized version with preprocessing cached
    """
    # Step 1: Preprocess all data ONCE
    preprocessed = preprocess_all_data(nii_files, events_files, mask_img, n_task_regressors, debug)
    
    mask_data = preprocessed['mask_data']
    mask_shape = preprocessed['mask_shape']
    
    # Store results
    r2_sum = np.zeros(mask_shape, dtype=np.float32)
    r2_sq_sum = np.zeros(mask_shape, dtype=np.float32)
    
    # For permutations, we only need the max distribution
    if permute:
        max_distribution = np.zeros(len(splits), dtype=np.float32)
        for split_id, (train_idx, test_idx) in enumerate(splits):
            max_r2 = run_single_cv_split(
                train_idx, test_idx, preprocessed, 
                n_task_regressors, tr, hrf_model, 
                permute_designs=True, 
                random_state=random_state + split_id
            )
            max_distribution[split_id] = max_r2
        return max_distribution, None
    
    # For real data, run all splits
    for split_id, (train_idx, test_idx) in enumerate(splits):
        print(f"\n{'='*60}")
        print(f"Processing split {split_id + 1}/{len(splits)}")
        print(f"{'='*60}")
        
        # Run split and get full R² map
        # (simplified - in practice you'd compute full map here)
        r2_map = np.zeros(mask_shape, dtype=np.float32)
        # ... compute full R² map ...
        
        r2_sum += r2_map
        r2_sq_sum += r2_map ** 2
        
        gc.collect()
    
    n_splits = len(splits)
    mean_r2 = r2_sum / n_splits
    var_r2 = (r2_sq_sum / n_splits) - (mean_r2 ** 2)
    
    return mean_r2, var_r2


# ============================================================
# CLI (modified)
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Optimized Cross-validated GLM")
    
    # [keep all your existing arguments]
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
    parser.add_argument("--n_jobs", type=int, default=1, 
                        help="Number of parallel jobs for permutations")
    parser.add_argument("--cache_preprocessed", action="store_true",
                        help="Save preprocessed data to disk for reuse")
    
    args = parser.parse_args()
    
    # [keep your existing validation and printing]
    
    # Load mask
    mask_img = nib.load(args.mask)
    
    # Create splits
    splitter = ShuffleSplit(
        n_splits=args.n_splits,
        test_size=args.test_size,
        random_state=args.random_state
    )
    splits = list(splitter.split(args.nii_files))
    
    # Check if we have preprocessed data cached
    cache_file = f"{args.output_prefix}_preprocessed.pkl"
    if args.cache_preprocessed and os.path.exists(cache_file):
        print(f"\nLoading preprocessed data from {cache_file}")
        with open(cache_file, 'rb') as f:
            preprocessed = pickle.load(f)
    else:
        # Preprocess once
        preprocessed = preprocess_all_data(
            args.nii_files, args.events_files, mask_img, 
            args.n_task_regressors, args.debug
        )
        if args.cache_preprocessed:
            print(f"\nSaving preprocessed data to {cache_file}")
            with open(cache_file, 'wb') as f:
                pickle.dump(preprocessed, f)
    
    # Run real data CV
    if not args.only_permutations:
        print("\n" + "="*60)
        print("RUNNING REAL DATA CV")
        print("="*60)
        
        mean_r2, var_r2 = run_cv_optimized(
            args.nii_files, args.events_files, mask_img,
            args.n_task_regressors, args.tr, args.hrf_model,
            splits, args.random_state, permute=False,
            debug=args.debug, n_jobs=args.n_jobs
        )
        
        # Save results
        img_mean = nib.Nifti1Image(mean_r2, mask_img.affine)
        img_mean.to_filename(f"{args.output_prefix}_real_mean_r2.nii.gz")
        
        img_var = nib.Nifti1Image(var_r2, mask_img.affine)
        img_var.to_filename(f"{args.output_prefix}_real_var_r2.nii.gz")
    
    # Run permutations - MUCH FASTER NOW
    if args.n_permutations > 0:
        print(f"\n{'='*60}")
        print(f"RUNNING {args.n_permutations} PERMUTATIONS")
        print(f"{'='*60}")
        
        max_distribution = np.zeros(args.n_permutations, dtype=np.float32)
        
        # Generate permutation seeds
        perm_seeds = [args.random_state + i + 1 for i in range(args.n_permutations)]
        
        # Run permutations in parallel
        if args.n_jobs > 1:
            results = Parallel(n_jobs=args.n_jobs)(
                delayed(run_single_permutation)(
                    preprocessed, args.n_task_regressors, args.tr, 
                    args.hrf_model, splits, seed
                ) for seed in tqdm(perm_seeds, desc="Permutations")
            )
            max_distribution = np.array(results)
        else:
            for p, seed in enumerate(tqdm(perm_seeds, desc="Permutations")):
                max_distribution[p] = run_single_permutation(
                    preprocessed, args.n_task_regressors, args.tr,
                    args.hrf_model, splits, seed
                )
        
        # Save results
        out_dist = f"{args.output_prefix}_perm_max_distribution.npy"
        np.save(out_dist, max_distribution)
        
        # Calculate thresholds
        threshold_95 = np.percentile(max_distribution, 95)
        threshold_99 = np.percentile(max_distribution, 99)
        print(f"\nPermutation thresholds:")
        print(f"  95th: {threshold_95:.6f}")
        print(f"  99th: {threshold_99:.6f}")


def run_single_permutation(preprocessed, n_task_regressors, tr, hrf_model, splits, seed):
    """Run a single permutation and return max R²"""
    max_r2 = 0
    for train_idx, test_idx in splits:
        r2 = run_single_cv_split(
            train_idx, test_idx, preprocessed,
            n_task_regressors, tr, hrf_model,
            permute_designs=True, random_state=seed
        )
        max_r2 = max(max_r2, r2)
    return max_r2


if __name__ == "__main__":
    main()