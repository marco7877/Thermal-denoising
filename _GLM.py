#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.glm.first_level import FirstLevelModel, make_first_level_design_matrix
from nilearn.masking import apply_mask, unmask
from nilearn.image import concat_imgs, get_data
from joblib import Parallel, delayed
from tqdm.auto import tqdm


# ============================================================
# Split generator (explicitly controlled by n_splits)
# ============================================================

def generate_splits(n_runs, n_splits, n_hold_out, random_state):

    rng = np.random.RandomState(random_state)
    indices = np.arange(n_runs)

    splits = []

    for _ in range(n_splits):

        test_idx = rng.choice(
            indices,
            size=(n_runs // n_hold_out),
            replace=False
        )

        train_idx = np.setdiff1d(indices, test_idx)

        splits.append((train_idx, test_idx))

    return splits


# ============================================================
# Core CV computation
# ============================================================

def run_cv(
    nii_files,
    events_files,
    mask_img,
    TR,
    HRF,
    n_regressors,
    splits,
    random_state,
    n_jobs,
    permute=False
):

    rng = np.random.RandomState(random_state)
    n_runs = len(nii_files)

    masked_data = [apply_mask(f, mask_img) for f in nii_files]
    n_timepoints = masked_data[0].shape[0]
    frame_times = np.arange(n_timepoints) * TR

    design_matrices = []

    for i in range(n_runs):
        tmp = pd.read_csv(events_files[i], sep="\t")
        events = tmp.loc[tmp["trial_type"] != "baseline"]

        dm = make_first_level_design_matrix(
            frame_times,
            events,
            drift_model="polynomial",
            drift_order=4,
            hrf_model=HRF
        ).fillna(0)

        design_matrices.append(dm)

    def run_split(train, test):

        train_imgs = [nii_files[i] for i in train]
        test_data = np.concatenate([masked_data[i] for i in test], axis=0)

        design_train = pd.concat(
            [design_matrices[i] for i in train],
            ignore_index=True
        )

        design_test = pd.concat(
            [design_matrices[i] for i in test],
            ignore_index=True
        )

        if permute:
            perm_idx = rng.permutation(len(design_train))
            design_train.iloc[:, :n_regressors] = \
                design_train.iloc[perm_idx, :n_regressors].values

        glm = FirstLevelModel(
            t_r=TR,
            mask_img=mask_img,
            standardize=False,
            signal_scaling=False,
            hrf_model=HRF,
            minimize_memory=True
        )

        glm = glm.fit(concat_imgs(train_imgs), design_matrices=design_train)

        contrast = np.eye(n_regressors, len(design_train.columns))

        betas = get_data(
            glm.compute_contrast(
                contrast,
                output_type="effect_size"
            )
        )

        test_design = design_test.values[:, :n_regressors]

        predicted = np.tensordot(
            betas,
            test_design,
            axes=([3], [1])
        )

        predicted = predicted.reshape(-1, predicted.shape[-1])
        y = test_data.T

        # Correct Pearson R²
        x = predicted - predicted.mean(axis=1, keepdims=True)
        y = y - y.mean(axis=1, keepdims=True)

        num = np.sum(x * y, axis=1)
        den = np.sqrt(np.sum(x**2, axis=1) * np.sum(y**2, axis=1))

        r = num / (den + 1e-8)

        return r**2

    # Parallel execution with progress bar
    with Parallel(n_jobs=n_jobs) as parallel:
        results = parallel(
            delayed(run_split)(train, test)
            for train, test in tqdm(splits, desc="CV splits")
        )

    results = np.array(results)
    return results.mean(axis=0)


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Cross-validated GLM with controlled splits and permutations"
    )

    parser.add_argument("--nii_files", nargs="+", required=True)
    parser.add_argument("--events_files", nargs="+", required=True)
    parser.add_argument("--mask", required=True)
    parser.add_argument("--output_prefix", required=True)

    parser.add_argument("--TR", type=float, default=2.0)
    parser.add_argument("--HRF", default="spm")
    parser.add_argument("--n_regressors", type=int, default=6)

    # 🔥 Explicit split control
    parser.add_argument("--n_splits", type=int, required=True)
    parser.add_argument("--n_hold_out", type=int, default=3)

    parser.add_argument("--n_permutations", type=int, default=0)
    parser.add_argument("--random_state", type=int, default=42)
    parser.add_argument("--n_jobs", type=int, default=-1)

    args = parser.parse_args()

    n_runs = len(args.nii_files)

    splits = generate_splits(
        n_runs=n_runs,
        n_splits=args.n_splits,
        n_hold_out=args.n_hold_out,
        random_state=args.random_state
    )

    mask_img = nib.load(args.mask)

    print("Running real-data CV...")
    real_r2 = run_cv(
        args.nii_files,
        args.events_files,
        mask_img,
        args.TR,
        args.HRF,
        args.n_regressors,
        splits,
        args.random_state,
        args.n_jobs,
        permute=False
    )

    unmask(real_r2, mask_img).to_filename(
        f"{args.output_prefix}_real_mean_r2.nii.gz"
    )

    # --------------------------------------------------------
    # Permutations
    # --------------------------------------------------------

    if args.n_permutations > 0:

        print(f"Running {args.n_permutations} permutations...")

        perm_maps = []

        for p in tqdm(range(args.n_permutations), desc="Permutations"):

            perm_r2 = run_cv(
                args.nii_files,
                args.events_files,
                mask_img,
                args.TR,
                args.HRF,
                args.n_regressors,
                splits,
                args.random_state + p + 1,
                args.n_jobs,
                permute=True
            )

            perm_maps.append(perm_r2)

        perm_maps = np.array(perm_maps)

        unmask(perm_maps.mean(axis=0), mask_img).to_filename(
            f"{args.output_prefix}_perm_mean_r2.nii.gz"
        )

        np.save(
            f"{args.output_prefix}_perm_max_distribution.npy",
            perm_maps.max(axis=1)
        )

    print("Finished.")


if __name__ == "__main__":
    main()
