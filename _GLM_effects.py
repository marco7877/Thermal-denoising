#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from nibabel import load, Nifti1Image
from nilearn.image import concat_imgs, get_data, new_img_like, resample_to_img
from nilearn.glm.first_level import FirstLevelModel
import warnings
warnings.filterwarnings('ignore')

# Optional import for plotting
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_and_preprocess_runs(nii_files, debug=False):
    """Load runs, resample to first run if needed (no scaling)."""
    all_runs_imgs = []
    run_timepoints = []
    reference_img = None

    for f in nii_files:
        img = load(f)
        if reference_img is None:
            reference_img = img
        else:
            if not np.allclose(img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
                img = resample_to_img(img, reference_img, interpolation='continuous')
        run_timepoints.append(img.shape[-1])
        all_runs_imgs.append(img)

    return all_runs_imgs, reference_img, run_timepoints


def load_design_matrices(events_files, run_timepoints, n_task_regressors):
    """Load design matrices and verify timepoint match."""
    designs = []
    for i, f in enumerate(events_files):
        df = pd.read_csv(f)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        if len(df) != run_timepoints[i]:
            raise ValueError(f"Timepoint mismatch for run {i}: fMRI={run_timepoints[i]}, Design={len(df)}")
        designs.append(df)
    return designs


def build_block_design(designs, n_task_regressors):
    """Build block‑diagonal design matrix with run‑specific nuisance columns."""
    blocks = []
    for run_idx, df in enumerate(designs):
        task_part = df.iloc[:, :n_task_regressors]
        nuisance_part = df.iloc[:, n_task_regressors:]
        rename_dict = {col: f"{col}_run{run_idx}" for col in nuisance_part.columns}
        nuisance_part = nuisance_part.rename(columns=rename_dict)
        run_design = pd.concat([task_part, nuisance_part], axis=1)
        blocks.append(run_design)
    return pd.concat(blocks, axis=0, ignore_index=True).fillna(0)


def save_design_matrix_image(design_matrix, output_file):
    """Save a heatmap of the design matrix to a PNG file."""
    if not HAS_MPL:
        print("Warning: matplotlib not installed. Cannot save design matrix image.")
        return

    plt.figure(figsize=(12, 8))
    plt.imshow(design_matrix.values, aspect='auto', cmap='RdBu_r', interpolation='none')
    plt.colorbar(label='Regressor value')
    plt.title('Design matrix')
    plt.xlabel('Regressors')
    plt.ylabel('Timepoints')
    n_cols = design_matrix.shape[1]
    if n_cols > 50:
        plt.xticks(np.arange(0, n_cols, step=max(1, n_cols//20)), rotation=90, fontsize=8)
    else:
        plt.xticks(np.arange(n_cols), design_matrix.columns, rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Design matrix image saved to {output_file}")


def check_mask(mask_img, reference_img=None):
    """Check that mask has positive voxels and optionally resample to reference."""
    mask_data = mask_img.get_fdata()
    if np.sum(mask_data) == 0:
        raise ValueError("Mask has zero voxels. Check your mask file.")
    
    if reference_img is not None:
        # Check if affine matches; if not, warn and optionally resample
        if not np.allclose(mask_img.affine, reference_img.affine, rtol=1e-3, atol=1e-3):
            print("Warning: Mask affine differs from reference image. Resampling mask to reference space.")
            mask_img = resample_to_img(mask_img, reference_img, interpolation='nearest')
            # After resampling, check again
            if np.sum(mask_img.get_fdata()) == 0:
                raise ValueError("After resampling, mask has zero voxels. Check mask and reference alignment.")
    return mask_img


def apply_mask_to_image(img, mask_img, fill_value=0):
    """Set voxels outside the mask to fill_value (default 0)."""
    data = get_data(img)
    mask_data = mask_img.get_fdata().astype(bool)
    # Ensure mask has same spatial dimensions as data
    if data.shape[:-1] != mask_data.shape:
        raise ValueError(f"Mask shape {mask_data.shape} does not match image spatial shape {data.shape[:-1]}")
    data_out = np.where(mask_data[..., np.newaxis], data, fill_value)
    return new_img_like(img, data_out)


def main():
    parser = argparse.ArgumentParser(description="Simple first‑level GLM")
    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--n_task_regressors", type=int, required=True)
    parser.add_argument("--contrast", type=str,
                        default="face1-(bodylimb1+ES_SC1+ES_RW1+ES_FF1+ES_CS1)/5")
    parser.add_argument("--tr", type=float, default=2.0)
    parser.add_argument("--hrf_model", type=str, default="glover",
                        help="HRF model: 'spm', 'glover', 'spm + derivative', etc.")
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--mask_output", action="store_true",
                        help="Apply mask to output images (set outside mask to 0)")
    parser.add_argument("--fill_nan", action="store_true",
                        help="Replace NaN values with 0 in output images")
    parser.add_argument("--save_design_matrix", action="store_true",
                        help="Save design matrix as PNG image (requires matplotlib)")
    args = parser.parse_args()

    # Validate inputs
    if len(args.nii_files) != len(args.events_files):
        raise ValueError("Number of NIfTI files and event files must match")

    # Load and preprocess runs (no scaling)
    all_runs_imgs, reference_img, run_timepoints = load_and_preprocess_runs(
        args.nii_files, debug=args.debug
    )

    # Load mask, optionally resample to first run's space
    mask_img = load(args.mask)
    mask_img = check_mask(mask_img, reference_img=reference_img)

    # Load design matrices
    designs = load_design_matrices(args.events_files, run_timepoints, args.n_task_regressors)

    # Build block‑diagonal design
    design_matrix = build_block_design(designs, args.n_task_regressors)

    # Print design matrix details
    print("\n" + "="*60)
    print("DESIGN MATRIX")
    print("="*60)
    print(f"Shape: {design_matrix.shape} (rows = total timepoints, columns = regressors)")
    print(f"Columns: {design_matrix.columns.tolist()}")
    print(f"First 5 rows:\n{design_matrix.head()}")
    print("="*60 + "\n")

    # Check design matrix rank
    rank = np.linalg.matrix_rank(design_matrix.values)
    print(f"Design matrix rank: {rank} (full rank would be {design_matrix.shape[1]})")
    if rank < design_matrix.shape[1]:
        print("WARNING: Design matrix is not full rank. This may cause estimation problems.")
        # Optionally, you could drop collinear columns here.

    # Save design matrix image if requested
    if args.save_design_matrix:
        design_img_file = f"{args.output_prefix}_design_matrix.png"
        save_design_matrix_image(design_matrix, design_img_file)

    # Concatenate runs
    concat_img = concat_imgs(all_runs_imgs)

    # Optional: check variance of masked data to identify constant voxels
    mask_data = mask_img.get_fdata().astype(bool)
    concat_data = get_data(concat_img)
    if mask_data.sum() > 0:
        masked_data = concat_data[mask_data, :]
        variance = np.var(masked_data, axis=1)
        n_const = np.sum(variance < 1e-6)
        print(f"Number of constant voxels within mask: {n_const} / {mask_data.sum()} ({100*n_const/mask_data.sum():.2f}%)")

    # Fit GLM
    print("Fitting GLM...")
    glm = FirstLevelModel(
        t_r=args.tr,
        mask_img=None,          # use whole image
        standardize=True,
        signal_scaling=False,
        hrf_model=args.hrf_model,
        minimize_memory=True,
    )
    glm = glm.fit(concat_img, design_matrices=design_matrix)

    # Compute contrast
    print(f"Computing contrast: {args.contrast}")
    t_img = glm.compute_contrast(args.contrast, stat_type="t")
    beta_img = glm.compute_contrast(args.contrast, stat_type="effect_size")

    # Optionally apply mask and/or fill NaN
    if args.mask_output:
        t_img = apply_mask_to_image(t_img, mask_img)
        beta_img = apply_mask_to_image(beta_img, mask_img)

    if args.fill_nan:
        t_data = get_data(t_img)
        t_data = np.nan_to_num(t_data, nan=0.0)
        t_img = new_img_like(t_img, t_data)
        beta_data = get_data(beta_img)
        beta_data = np.nan_to_num(beta_data, nan=0.0)
        beta_img = new_img_like(beta_img, beta_data)

    # Diagnostic print
    t_data = get_data(t_img)
    beta_data = get_data(beta_img)
    mask_data = mask_img.get_fdata().astype(bool)

    t_masked = t_data[mask_data]
    beta_masked = beta_data[mask_data]

    n_nan_t = np.isnan(t_masked).sum()
    n_nan_beta = np.isnan(beta_masked).sum()
    n_voxels = mask_data.sum()

    print(f"\nDiagnostics within mask ({n_voxels} voxels):")
    print(f"  t‑statistic: NaNs = {n_nan_t} ({100*n_nan_t/n_voxels:.2f}%)")
    if not np.all(np.isnan(t_masked)):
        print(f"    min/max = {np.nanmin(t_masked):.4f} / {np.nanmax(t_masked):.4f}")
        print(f"    mean ± std = {np.nanmean(t_masked):.4f} ± {np.nanstd(t_masked):.4f}")
    print(f"  beta (effect): NaNs = {n_nan_beta} ({100*n_nan_beta/n_voxels:.2f}%)")
    if not np.all(np.isnan(beta_masked)):
        print(f"    min/max = {np.nanmin(beta_masked):.4f} / {np.nanmax(beta_masked):.4f}")
        print(f"    mean ± std = {np.nanmean(beta_masked):.4f} ± {np.nanstd(beta_masked):.4f}")

    # Save outputs
    t_out = f"{args.output_prefix}_t.nii.gz"
    beta_out = f"{args.output_prefix}_beta.nii.gz"
    t_img.to_filename(t_out)
    beta_img.to_filename(beta_out)
    print(f"\nSaved t‑statistic map: {t_out}")
    print(f"Saved beta map: {beta_out}")

    print("Done.")


if __name__ == "__main__":
    main()
