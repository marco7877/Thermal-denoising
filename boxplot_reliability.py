#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI tool for generating split violin plots comparing two pipelines across base methods and subjects.
Uses seaborn's split=True to show left/right halves for the two pipelines.
"""

import numpy as np
import pandas as pd
import os
import argparse
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns
from nilearn.image import load_img
from nilearn.masking import apply_mask

#####################################################################################
###### Functions ####################################################################
#####################################################################################

def find_mask_for_subject(directory, mask_pattern, subject):
    """Find mask file for a specific subject using pattern."""
    mask_file_pattern = mask_pattern.replace("{subject}", subject)
    full_pattern = os.path.join(directory, mask_file_pattern)
    matching_files = glob.glob(full_pattern)
    
    if matching_files:
        print(f"  Found mask for {subject}: {os.path.basename(matching_files[0])}")
        return matching_files[0]
    return None


def load_and_reshape_nifti(file_path, mask_img=None):
    """Load a nifti file and reshape to 1D array, optionally applying a mask."""
    print(f"    Loading: {os.path.basename(file_path)}")
    img = load_img(file_path)
    
    if mask_img is not None:
        data = apply_mask(img, mask_img)
        if data.ndim > 1:
            data = data.flatten()
    else:
        data = np.array(img.dataobj).flatten()
    
    return data


def find_files_for_subject(directory, pattern, subject, full_methods):
    """Find files for a subject for all full method names (including pipeline suffix)."""
    file_dict = {}
    for method in full_methods:
        file_pattern = pattern.replace("{subject}", subject).replace("{method}", method)
        full_pattern = os.path.join(directory, file_pattern)
        matching_files = glob.glob(full_pattern)
        
        if matching_files:
            file_dict[method] = matching_files[0]
            print(f"  Found {method}: {os.path.basename(matching_files[0])}")
        else:
            print(f"  Warning: No file for {subject}, method {method}")
            file_dict[method] = None
    
    return file_dict


def collect_all_data(directory, pattern, subjects, base_methods,
                     pipeline1_suffix, pipeline2_suffix,
                     drop_zeros, mask_dir=None, mask_pattern=None):
    """
    Load all data, applying masks if provided.
    Returns a DataFrame with columns: subject, base_method, pipeline, tsnr
    """
    using_masks = mask_dir and mask_pattern
    all_data = []
    
    # Full method names for both pipelines
    pipeline1_methods = [base + pipeline1_suffix for base in base_methods]
    pipeline2_methods = [base + pipeline2_suffix for base in base_methods]
    all_full_methods = pipeline1_methods + pipeline2_methods
    
    for subject in subjects:
        print(f"\nProcessing {subject}")
        
        # Load mask if provided
        mask_img = None
        if using_masks:
            mask_file = find_mask_for_subject(mask_dir, mask_pattern, subject)
            if mask_file:
                mask_img = load_img(mask_file)
        
        # Find files
        tsnr_files = find_files_for_subject(directory, pattern, subject, all_full_methods)
        
        for full_method, file_path in tsnr_files.items():
            if file_path is None:
                continue
            
            data = load_and_reshape_nifti(file_path, mask_img)
            if drop_zeros:
                data = data[data != 0]
            
            # Determine base method and pipeline
            if full_method.endswith(pipeline1_suffix):
                base = full_method[:-len(pipeline1_suffix)]
                pipeline = pipeline1_suffix
            elif full_method.endswith(pipeline2_suffix):
                base = full_method[:-len(pipeline2_suffix)]
                pipeline = pipeline2_suffix
            else:
                continue  # should not happen
            
            df_temp = pd.DataFrame({
                'subject': subject,
                'base_method': base,
                'pipeline': pipeline,
                'tsnr': data
            })
            all_data.append(df_temp)
            print(f"    Added {len(data)} voxels for {base} ({pipeline})")
    
    if not all_data:
        print("ERROR: No data collected.")
        return None
    
    df_all = pd.concat(all_data, ignore_index=True)
    # Set categorical orders
    df_all['base_method'] = pd.Categorical(df_all['base_method'], categories=base_methods, ordered=True)
    df_all['pipeline'] = pd.Categorical(df_all['pipeline'], categories=[pipeline1_suffix, pipeline2_suffix], ordered=True)
    # Ensure subjects are ordered (use the order they appear in subjects list)
    df_all['subject'] = pd.Categorical(df_all['subject'], categories=subjects, ordered=True)
    
    return df_all


def compute_summary(df_all):
    """Compute mean, median, std per subject, base_method, pipeline."""
    summary = df_all.groupby(['subject', 'base_method', 'pipeline'], observed=False)['tsnr'] \
                    .agg(['mean', 'median', 'std']).reset_index()
    summary['ymin'] = summary['mean'] - summary['std']
    summary['ymax'] = summary['mean'] + summary['std']
    return summary


def plot_split_violin_across_subjects(df_all, summary, base_methods, pipeline1_suffix, pipeline2_suffix,
                                       output_path, ylim, show_points, alpha, dpi=150):
    """
    Create a split violin plot with:
        x-axis: subject_method combination (subject + method)
        split by pipeline (left/right)
    """
    # Create a combined column for x-axis: subject + method
    df_all['subject_method'] = df_all['subject'].astype(str) + "_" + df_all['base_method'].astype(str)
    # Order subject_method: first all methods for subject1 in base_methods order, then subject2, etc.
    ordered_subject_method = []
    for subj in df_all['subject'].cat.categories:
        for method in base_methods:
            ordered_subject_method.append(f"{subj}_{method}")
    df_all['subject_method'] = pd.Categorical(df_all['subject_method'], categories=ordered_subject_method, ordered=True)
    
    # Prepare summary with same column
    summary['subject_method'] = summary['subject'].astype(str) + "_" + summary['base_method'].astype(str)
    summary['subject_method'] = pd.Categorical(summary['subject_method'], categories=ordered_subject_method, ordered=True)
    
    # Set style
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Split violinplot
    sns.violinplot(data=df_all, x='subject_method', y='tsnr', hue='pipeline',
                   split=True, inner='quart', palette='Set2',
                   linewidth=1, ax=ax)
    
    # Overlay mean ± SD points
    # We need to compute x positions for each subject_method and pipeline combination
    # Since we have two pipelines, we need to offset the points.
    # seaborn internally uses 0,1 for x positions, but after plotting we can extract positions.
    # Simpler: use matplotlib to plot points with manual x positions.
    # We'll compute numeric positions based on the order.
    
    # Get unique x tick positions (one per subject_method)
    x_ticks = ordered_subject_method
    # Map each subject_method to its x index
    x_map = {sm: i for i, sm in enumerate(x_ticks)}
    
    # For each pipeline, we need to shift left/right a bit
    # In split violins, left is usually pipeline1, right pipeline2. We'll use dodge of 0.2.
    dodge = 0.2
    pipeline_pos = {pipeline1_suffix: -dodge, pipeline2_suffix: dodge}
    
    # Plot mean points
    for pipeline in [pipeline1_suffix, pipeline2_suffix]:
        sub = summary[summary['pipeline'] == pipeline]
        x_pos = [x_map[sm] + pipeline_pos[pipeline] for sm in sub['subject_method']]
        ax.scatter(x_pos, sub['mean'], color='black', s=40, zorder=5,
                   label=f'Mean ({pipeline})' if pipeline == pipeline1_suffix else "")
        # Error bars
        ax.errorbar(x_pos, sub['mean'], yerr=sub['std'], fmt='none',
                    ecolor='black', capsize=3, elinewidth=1, zorder=4)
    
    # Add individual points if requested (sample)
    if show_points:
        # Sample up to 1000 points per group
        sampled = df_all.groupby(['subject_method', 'pipeline'], observed=False).apply(
            lambda x: x.sample(min(1000, len(x)))
        ).reset_index(drop=True)
        # Convert x categories to numeric positions
        sampled['x_numeric'] = sampled['subject_method'].map(x_map)
        # Shift by pipeline
        sampled['x_numeric'] = sampled.apply(
            lambda row: row['x_numeric'] + pipeline_pos[row['pipeline']], axis=1)
        ax.scatter(sampled['x_numeric'], sampled['tsnr'], alpha=alpha,
                   s=2, c='gray', label='_nolegend_')
    
    # Beautify
    ax.set_ylim(ylim)
    ax.set_xlabel("Subject - Method")
    ax.set_ylabel("tSNR")
    ax.set_title(f"tSNR distribution: {pipeline1_suffix} (left) vs {pipeline2_suffix} (right)\nMean ± SD shown as black points & error bars")
    
    # Rotate x labels
    plt.xticks(rotation=45, ha='right')
    # Move legend to right
    ax.legend(loc='upper right', title='Pipeline')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    print(f"Saved across-subjects plot to {os.path.basename(output_path)}")


def plot_split_violin_per_subject(df_all, summary, base_methods, pipeline1_suffix, pipeline2_suffix,
                                  output_dir, ylim, show_points, alpha, dpi=150, overwrite=True):
    """
    For each subject, create a split violin plot with:
        x-axis: base_method
        split by pipeline (left/right)
    """
    results = []
    subjects = df_all['subject'].cat.categories
    for subject in subjects:
        output_path = os.path.join(output_dir, f"{subject}_split_violin.png")
        if os.path.exists(output_path) and not overwrite:
            print(f"File exists. Skipping {subject}.")
            results.append(output_path)
            continue
        
        df_sub = df_all[df_all['subject'] == subject].copy()
        if df_sub.empty:
            continue
        
        # Ensure base_method order
        df_sub['base_method'] = pd.Categorical(df_sub['base_method'], categories=base_methods, ordered=True)
        # Summary for this subject
        sub_summary = summary[summary['subject'] == subject].copy()
        sub_summary['base_method'] = pd.Categorical(sub_summary['base_method'], categories=base_methods, ordered=True)
        
        sns.set_style("whitegrid")
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Split violinplot
        sns.violinplot(data=df_sub, x='base_method', y='tsnr', hue='pipeline',
                       split=True, inner='quart', palette='Set2',
                       linewidth=1, ax=ax)
        
        # Overlay mean ± SD points
        # Get numeric positions for each base_method (0,1,2,...)
        x_positions = np.arange(len(base_methods))
        dodge = 0.2
        pipeline_pos = {pipeline1_suffix: -dodge, pipeline2_suffix: dodge}
        
        for pipeline in [pipeline1_suffix, pipeline2_suffix]:
            subp = sub_summary[sub_summary['pipeline'] == pipeline]
            # Map base_method to its numeric index
            x_vals = [x_positions[base_methods.index(m)] + pipeline_pos[pipeline] for m in subp['base_method']]
            ax.scatter(x_vals, subp['mean'], color='black', s=40, zorder=5,
                       label=f'Mean ({pipeline})' if pipeline == pipeline1_suffix else "")
            ax.errorbar(x_vals, subp['mean'], yerr=subp['std'], fmt='none',
                        ecolor='black', capsize=3, elinewidth=1, zorder=4)
        
        # Add individual points if requested
        if show_points:
            sampled = df_sub.groupby(['base_method', 'pipeline'], observed=False).apply(
                lambda x: x.sample(min(1000, len(x)))
            ).reset_index(drop=True)
            # Map base_method to index and shift
            sampled['x_numeric'] = sampled['base_method'].map({m: i for i, m in enumerate(base_methods)})
            sampled['x_numeric'] = sampled.apply(
                lambda row: row['x_numeric'] + pipeline_pos[row['pipeline']], axis=1)
            ax.scatter(sampled['x_numeric'], sampled['tsnr'], alpha=alpha,
                       s=2, c='gray', label='_nolegend_')
        
        ax.set_ylim(ylim)
        ax.set_xlabel("Method")
        ax.set_ylabel("tSNR")
        ax.set_title(f"{subject}: tSNR distribution\n{pipeline1_suffix} (left) vs {pipeline2_suffix} (right)")
        plt.xticks(rotation=45, ha='right')
        ax.legend(loc='upper right', title='Pipeline')
        plt.tight_layout()
        plt.savefig(output_path, dpi=dpi)
        plt.close()
        print(f"Saved {os.path.basename(output_path)}")
        results.append(output_path)
    
    return results


def find_available_subjects(directory, pattern, base_methods, pipeline1_suffix, pipeline2_suffix):
    """Find subjects that have all required files (all base_methods with both pipelines)."""
    full_methods = [base + pipeline1_suffix for base in base_methods] + [base + pipeline2_suffix for base in base_methods]
    search_pattern = pattern.replace("{subject}", "*").replace("{method}", "*")
    all_files = glob.glob(os.path.join(directory, search_pattern))
    
    subjects = set()
    for f in all_files:
        match = re.search(r'(sub-[^_]+)', os.path.basename(f))
        if match:
            subjects.add(match.group(1))
    
    valid = []
    for s in sorted(subjects):
        if all(glob.glob(os.path.join(directory,
               pattern.replace("{subject}", s).replace("{method}", m))) 
               for m in full_methods):
            valid.append(s)
    
    return valid


#####################################################################################
###### Main #########################################################################
#####################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Generate split violin plots comparing two pipelines across base methods and subjects. "
                    "Uses seaborn's split=True.")
    
    parser.add_argument("--directory", type=str, required=True,
                       help="Directory with nifti files")
    
    parser.add_argument("--pattern", type=str, required=True,
                       help="File pattern with {subject} and {method}. Example: '{subject}_task_{method}.nii.gz'")
    
    parser.add_argument("--base_methods", type=str, nargs="+", required=True,
                       help="Base method names (e.g., method1 method2 method3) – without pipeline suffix")
    
    parser.add_argument("--pipeline1_suffix", type=str, required=True,
                       help="Suffix for first pipeline (e.g., '01') – appended to base methods")
    
    parser.add_argument("--pipeline2_suffix", type=str, required=True,
                       help="Suffix for second pipeline (e.g., '101') – appended to base methods")
    
    parser.add_argument("--mask_dir", type=str, default=None,
                       help="Directory with subject masks")
    
    parser.add_argument("--mask_pattern", type=str, default=None,
                       help="Mask pattern with {subject}. Example: '{subject}_mask.nii.gz'")
    
    # Subject selection
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--subjects", type=str, nargs="+",
                      help="Subject IDs (e.g., sub-001 sub-002)")
    group.add_argument("--auto_subjects", action="store_true",
                      help="Auto-detect subjects that have all required files")
    
    parser.add_argument("--output_dir", type=str, default=None,
                       help="Output directory (default: same as input)")
    
    parser.add_argument("--ylim", type=float, nargs=2, default=[0, 400],
                       help="Y-axis limits. Default: 0 400")
    
    parser.add_argument("--alpha", type=float, default=0.3,
                       help="Point transparency. Default: 0.3")
    
    parser.add_argument("--no_drop_zeros", action="store_true",
                       help="Keep zero values")
    
    parser.add_argument("--show_points", action="store_true",
                       help="Show individual points (sampled)")
    
    parser.add_argument("--no_overwrite", action="store_true",
                       help="Don't overwrite existing files")
    
    args = parser.parse_args()
    
    # Validate
    if (args.mask_dir and not args.mask_pattern) or (args.mask_pattern and not args.mask_dir):
        parser.error("Both --mask_dir and --mask_pattern are required together")
    
    if len(args.base_methods) < 1:
        parser.error("At least one base method required")
    
    # Determine subjects
    if args.auto_subjects:
        print(f"\nAuto-detecting subjects...")
        subjects = find_available_subjects(args.directory, args.pattern,
                                           args.base_methods,
                                           args.pipeline1_suffix, args.pipeline2_suffix)
        if not subjects:
            print("ERROR: No subjects found with all required files.")
            return
        print(f"Found: {', '.join(subjects)}")
    else:
        subjects = args.subjects
        print(f"\nSubjects: {', '.join(subjects)}")
    
    # Set output directory
    output_dir = args.output_dir if args.output_dir else args.directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect data
    df_all = collect_all_data(
        directory=args.directory,
        pattern=args.pattern,
        subjects=subjects,
        base_methods=args.base_methods,
        pipeline1_suffix=args.pipeline1_suffix,
        pipeline2_suffix=args.pipeline2_suffix,
        drop_zeros=not args.no_drop_zeros,
        mask_dir=args.mask_dir,
        mask_pattern=args.mask_pattern
    )
    
    if df_all is None:
        return
    
    # Compute summary statistics
    summary = compute_summary(df_all)
    summary_path = os.path.join(output_dir, "tsnr_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"\nSaved summary to {os.path.basename(summary_path)}")
    
    # Across-subjects plot
    across_path = os.path.join(output_dir, "split_violin_across_subjects.png")
    if not (os.path.exists(across_path) and args.no_overwrite):
        plot_split_violin_across_subjects(
            df_all, summary,
            base_methods=args.base_methods,
            pipeline1_suffix=args.pipeline1_suffix,
            pipeline2_suffix=args.pipeline2_suffix,
            output_path=across_path,
            ylim=tuple(args.ylim),
            show_points=args.show_points,
            alpha=args.alpha
        )
    else:
        print(f"Across-subjects plot exists. Skipping.")
    
    # Per-subject plots
    per_subject_results = plot_split_violin_per_subject(
        df_all, summary,
        base_methods=args.base_methods,
        pipeline1_suffix=args.pipeline1_suffix,
        pipeline2_suffix=args.pipeline2_suffix,
        output_dir=output_dir,
        ylim=tuple(args.ylim),
        show_points=args.show_points,
        alpha=args.alpha,
        overwrite=not args.no_overwrite
    )
    
    print(f"\n{'='*60}")
    print(f"DONE! Generated {len(per_subject_results)} per-subject plots and 1 across-subjects plot.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
