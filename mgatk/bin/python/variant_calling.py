#!/usr/bin/python

###################################################
# Call mito variants
###################################################

import sys
import gzip
import glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Parse arguments
MGATK_OUT_DIR = sys.argv[1]
sample_prefix = sys.argv[2] 
mito_length = int(sys.argv[3])
low_coverage_threshold = int(sys.argv[4])
mito_genome = sys.argv[5]
EPSILON = 1e-10
letters = list('ATCG')


def load_mgatk_output(output_dir, mito_length):
    """Load and process base coverage files."""
    base_files = [glob.glob(output_dir + f'*.{nt}.txt.gz')[0] for nt in 'ATCG']
    base_coverage_dict = dict()
    
    for i, nt in enumerate('ATCG'):
        cur_base_data = pd.read_csv(gzip.open(base_files[i]), header=None)
        positions = range(1, mito_length + 1)
        
        # Process forward and reverse strands
        fwd_pivot = cur_base_data[[0, 1, 2]].pivot_table(index=1, columns=0)
        fwd_pivot.columns = [x[1] for x in fwd_pivot.columns.values]
        rev_pivot = cur_base_data[[0, 1, 3]].pivot_table(index=1, columns=0)
        rev_pivot.columns = [x[1] for x in rev_pivot.columns.values]
        
        # Create complete DataFrames
        fwd_base_df = pd.DataFrame(0, index=fwd_pivot.index, columns=positions)
        rev_base_df = pd.DataFrame(0, index=rev_pivot.index, columns=positions)
        fwd_base_df.update(fwd_pivot)
        rev_base_df.update(rev_pivot)
        
        base_coverage_dict[nt] = (fwd_base_df, rev_base_df)
    
    return base_coverage_dict


def gather_possible_variants(base_coverage_dict, reference_file):
    """Identify potential variants from coverage data."""
    # Sum across cells and strands for each base and position
    aggregated_genotype = pd.DataFrame(
        np.zeros((4, mito_length)), 
        index=list('ATCG'),
        columns=np.arange(1, mito_length + 1)
    )
    
    for nt in base_coverage_dict:
        fwd_base_df, rev_base_df = base_coverage_dict[nt]
        fwd_base_sum, rev_base_sum = fwd_base_df.sum(), rev_base_df.sum()
        
        # Ignore sequencing artifacts
        masking = ~((fwd_base_sum > 0) & (rev_base_sum > 0))
        fwd_base_sum[masking], rev_base_sum[masking] = 0, 0
        
        # Sum across strands
        aggregated_genotype.loc[nt, :] = fwd_base_sum + rev_base_sum
    
    # Create reference set
    ref_set = [x.strip().split() for x in open(reference_file, 'r').readlines()]
    ref_N_positions = [int(x[0]) for x in ref_set if x[1].upper() not in letters]
    ref_set = set([(int(x[0]), x[1].upper()) for x in ref_set if x[1].upper() in letters])
    ref_dict = dict(ref_set)
    
    # Create observed set
    non_zero_idx = np.where(aggregated_genotype > 0)
    non_zero_bases = [letters[i] for i in non_zero_idx[0]]
    non_zero_pos = [int(i + 1) for i in non_zero_idx[1]]
    observed_set = set([
        (pos, base) for pos, base in zip(non_zero_pos, non_zero_bases) 
        if pos not in ref_N_positions
    ])
    
    # Determine variants
    variant_set = observed_set - ref_set
    variants = sorted([
        (x[0], ref_dict[x[0]], x[1]) for x in list(variant_set)
    ], key=lambda x: x[0])
    
    return variants


# Load base coverage data
base_coverage_dict = load_mgatk_output(MGATK_OUT_DIR, mito_length)
cell_barcodes = base_coverage_dict['A'][0].index

# Calculate total coverage per position per cell
total_coverage = pd.DataFrame(
    np.zeros((len(cell_barcodes), mito_length)),
    index=cell_barcodes,
    columns=np.arange(1, mito_length + 1)
)

for nt in base_coverage_dict:
    total_coverage += base_coverage_dict[nt][0]
    total_coverage += base_coverage_dict[nt][1]

# Exclude low coverage cells
cell_barcodes = total_coverage.index[total_coverage.mean(axis=1) > low_coverage_threshold]
for nt in base_coverage_dict:
    base_coverage_dict[nt] = (
        base_coverage_dict[nt][0].loc[cell_barcodes, :],
        base_coverage_dict[nt][1].loc[cell_barcodes, :]
    )
total_coverage = total_coverage.loc[cell_barcodes, :]

# Call potential variants
variants = gather_possible_variants(
    base_coverage_dict,
    MGATK_OUT_DIR + mito_genome + '_refAllele.txt'
)
variant_names = [f'{x[0]}{x[1]}>{x[2]}' for x in variants]

# Build cell by variant tables for each strand
total_coverage_variant_df = []
fwd_cell_variant_df, rev_cell_variant_df = [], []

for i, var in enumerate(variants):
    var_name = variant_names[i]
    pos, base = var[0], var[2]
    total_coverage_variant_df.append(total_coverage[pos])
    fwd_cell_variant_df.append(base_coverage_dict[base][0][pos].values)
    rev_cell_variant_df.append(base_coverage_dict[base][1][pos].values)

# Convert to DataFrames
total_coverage_variant_df = pd.DataFrame(
    np.array(total_coverage_variant_df).T,
    index=cell_barcodes,
    columns=variant_names
)
fwd_cell_variant_df = pd.DataFrame(
    np.array(fwd_cell_variant_df).T,
    index=cell_barcodes,
    columns=variant_names
)
rev_cell_variant_df = pd.DataFrame(
    np.array(rev_cell_variant_df).T,
    index=cell_barcodes,
    columns=variant_names
)
all_cell_variant_df = fwd_cell_variant_df + rev_cell_variant_df

# Calculate heteroplasmic ratio
heteroplasmic_df = all_cell_variant_df / (total_coverage_variant_df + EPSILON)

# Calculate strand correlation
mask_idx = (fwd_cell_variant_df + rev_cell_variant_df) == 0
fwd_cell_variant_df[mask_idx] = np.nan
rev_cell_variant_df[mask_idx] = np.nan
variant_strand_corr = fwd_cell_variant_df.corrwith(rev_cell_variant_df).round(3)

# Calculate VMR
variant_mean = all_cell_variant_df.sum() / (total_coverage_variant_df.sum() + EPSILON)
variant_var = heteroplasmic_df.var()
variant_vmr = variant_var / (variant_mean + EPSILON)

# Compute summary statistics
variant_stats = {
    'position': [x[0] for x in variants],
    'nucleotide': [f'{x[1]}>{x[2]}' for x in variants],
    'variant': variant_names,
    'vmr': variant_vmr,
    'mean': variant_mean,
    'variance': variant_var,
    'n_cells_conf_detected': ((fwd_cell_variant_df >= 2) & (rev_cell_variant_df >= 2)).sum(),
    'n_cells_over_5': (heteroplasmic_df >= 0.05).sum(),
    'n_cells_over_10': (heteroplasmic_df >= 0.1).sum(),
    'n_cells_over_20': (heteroplasmic_df >= 0.2).sum(),
    'n_cells_over_95': (heteroplasmic_df >= 0.95).sum(),
    'max_heteroplasmy': heteroplasmic_df.max(),
    'strand_correlation': variant_strand_corr,
    'mean_coverage': total_coverage_variant_df.mean()
}

# Create output DataFrame
variant_output = pd.DataFrame(variant_stats)

# Convert numeric columns to float64
numeric_cols = ['vmr', 'mean', 'variance', 'strand_correlation', 
                'mean_coverage', 'max_heteroplasmy']
variant_output[numeric_cols] = variant_output[numeric_cols].astype(np.float64)

# Filter variants
multi_cell_variants = variant_output[variant_output['n_cells_conf_detected'] >= 3]['variant']
heteroplasmic_df = heteroplasmic_df[multi_cell_variants]

# Generate plot
plt.figure(figsize=(10, 8))
plt.scatter(
    variant_output[variant_output['variant'].isin(multi_cell_variants)]['strand_correlation'],
    np.log10(variant_output[variant_output['variant'].isin(multi_cell_variants)]['vmr']),
    s=5
)
plt.axhline(np.log10(0.01), color='red', alpha=0.4, linestyle=':')
plt.axhline(np.log10(EPSILON), color='red', alpha=0.4, linestyle=':')
plt.axvline(0.65, color='red', alpha=0.4, linestyle=':')
plt.xlabel('strand correlation', fontsize=20)
plt.ylabel('log10(VMR)', fontsize=20)

# Save results
plt.savefig(MGATK_OUT_DIR + sample_prefix + '.vmr_strand_plot.png')
variant_output.to_csv(
    MGATK_OUT_DIR + sample_prefix + '.variant_stats.tsv.gz',
    sep='\t',
    compression='gzip',
    index=False
)
heteroplasmic_df.to_csv(
    MGATK_OUT_DIR + sample_prefix + '.cell_heteroplasmic_df.tsv.gz',
    sep='\t',
    compression='gzip'
)