# ==========================================================================================
# Script: Pairwise MPO and Docking Score Comparator with Benchmark Filtering
#
# GOAL:
# This script analyzes pairs of compounds (with identical titles) from an input CSV file 
# and calculates average and delta values for:
#   - MPO score (multi-parameter optimization score)
#   - Normalized docking score (binding affinity)
#
# It is designed to help evaluate compound pairs by:
#   - Ranking by average MPO and normalized docking score
#   - Comparing against a benchmark compound (e.g., cariprazine)
#   - Filtering and exporting only those compounds that exceed the benchmark on both metrics
#
# REQUIRED INPUT FILE:
# - A CSV file with at least the following columns:
#     - 'Title': Compound identifier (same for each enumerated pair)
#     - 'MPO_score': Multi-parameter optimization score (continuous 0–1)
#     - 'norm_docking_score': Normalized docking score (0–1)
#     - 'docking score': Raw docking score (e.g., GlideScore)
#
# USER INPUT:
# - Path to the CSV file
# - Benchmark compound name (case-insensitive partial match to any value in 'Title')
#
# OUTPUT FILES:
# 1. Sorted_by_Avg_MPO.csv
#     - All compounds sorted by average MPO score
# 2. Sorted_by_Avg_normDocking.csv
#     - All compounds sorted by average normalized docking score
#     - Includes MPO metrics following docking metrics for better visual comparison
# 3. Above_<benchmark>_Compounds.csv
#     - Only compounds with both Avg_MPO and Avg_norm_docking greater than the benchmark
#     - Appends the benchmark compound rows at the bottom for reference
#
# NOTES:
# - Rows with titles that do not appear exactly twice are passed through unchanged.
# - All sorting is in descending order for interpretability (higher = better).
# ==========================================================================================


import pandas as pd
import os

# Ask user for inputs
input_path = input("Enter the full path to the input CSV file: ").strip().replace("\\", "/").strip('"')
benchmark_name = input("Enter the benchmark compound name (e.g., cariprazine): ").strip().lower()

# Load dataset
df = pd.read_csv(input_path)

# Ensure required columns exist
required_columns = ['Title', 'MPO_score', 'norm_docking_score', 'docking score']
missing_cols = [col for col in required_columns if col not in df.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Group by Title and calculate averages and deltas
grouped = df.groupby('Title', sort=False)
processed_rows = []

for name, group in grouped:
    if len(group) != 2:
        processed_rows.append(group)
        continue
    row1, row2 = group.iloc[0].copy(), group.iloc[1].copy()
    avg_mpo = (row1['MPO_score'] + row2['MPO_score']) / 2
    delta_mpo = abs(row1['MPO_score'] - row2['MPO_score'])
    avg_dock = (row1['norm_docking_score'] + row2['norm_docking_score']) / 2
    delta_dock = abs(row1['norm_docking_score'] - row2['norm_docking_score'])
    for row in [row1, row2]:
        row['Avg_MPO'] = avg_mpo
        row['Delta_MPO'] = delta_mpo
        row['Avg_norm_docking'] = avg_dock
        row['Delta_norm_docking'] = delta_dock
    processed_rows.extend([row1, row2])

processed_df = pd.DataFrame(processed_rows)

# Sort and save by Avg MPO
cols = processed_df.columns.tolist()
for col in ['Avg_MPO', 'Delta_MPO', 'Avg_norm_docking', 'Delta_norm_docking', 'MPO_score', 'norm_docking_score']:
    if col in cols:
        cols.remove(col)
insert_idx = cols.index('Title') + 1
cols = (cols[:insert_idx] +
        ['Avg_MPO', 'Delta_MPO', 'MPO_score', 'Avg_norm_docking', 'Delta_norm_docking', 'norm_docking_score'] +
        cols[insert_idx:])
df_sorted_mpo = processed_df[cols].sort_values(by='Avg_MPO', ascending=False)
output_dir = os.path.dirname(input_path)
output_path_mpo = os.path.join(output_dir, "Sorted_by_Avg_MPO.csv")
df_sorted_mpo.to_csv(output_path_mpo, index=False)

# Sort and save by Avg norm docking
cols = processed_df.columns.tolist()
for col in ['Avg_MPO', 'Delta_MPO', 'MPO_score',
            'Avg_norm_docking', 'Delta_norm_docking',
            'norm_docking_score', 'docking score']:
    if col in cols:
        cols.remove(col)
insert_idx = cols.index('Title') + 1
cols = (cols[:insert_idx] +
        ['Avg_norm_docking', 'Delta_norm_docking', 'norm_docking_score', 'docking score',
         'Avg_MPO', 'Delta_MPO', 'MPO_score'] +
        cols[insert_idx:])
df_sorted_dock = processed_df[cols].sort_values(by='Avg_norm_docking', ascending=False)
output_path_dock = os.path.join(output_dir, "Sorted_by_Avg_normDocking.csv")
df_sorted_dock.to_csv(output_path_dock, index=False)

# Benchmark filtering
benchmark_row = processed_df[processed_df['Title'].str.lower().str.contains(benchmark_name, na=False)]
if benchmark_row.empty:
    print(f"\n⚠️ Benchmark compound '{benchmark_name}' not found.")
else:
    avg_mpo_benchmark = benchmark_row['Avg_MPO'].dropna().iloc[0]
    avg_norm_dock_benchmark = benchmark_row['Avg_norm_docking'].dropna().iloc[0]

    # Use Avg_norm_docking instead of norm_docking_score for comparison
    better_than_benchmark = processed_df[
        (processed_df['Avg_MPO'] > avg_mpo_benchmark) &
        (processed_df['Avg_norm_docking'] > avg_norm_dock_benchmark)
    ]
    
    final_output = pd.concat([better_than_benchmark, benchmark_row], ignore_index=True)
    # Reorder columns to place key fields at the top
    key_cols = ['Title', 'Avg_MPO', 'Avg_norm_docking', 'norm_docking_score', 'docking score']
    remaining_cols = [col for col in final_output.columns if col not in key_cols]
    final_output = final_output[key_cols + remaining_cols]
    output_path_filtered = os.path.join(output_dir, f"Above_{benchmark_name}_Compounds.csv")
    final_output.to_csv(output_path_filtered, index=False)

    # Final Summary
    print("\n✅ Processing complete.")
    unique_mpo = df_sorted_mpo[df_sorted_mpo['Avg_MPO'] > avg_mpo_benchmark]['Title'].nunique()
    unique_dock = df_sorted_dock[df_sorted_dock['Avg_norm_docking'] > avg_norm_dock_benchmark]['Title'].nunique()
    unique_both = better_than_benchmark['Title'].nunique()

    print(f"- File sorted by Avg MPO: {output_path_mpo}")
    print(f"  → {unique_mpo} unique compounds exceeded {benchmark_name.title()}'s Avg MPO")

    print(f"- File sorted by Avg norm docking score: {output_path_dock}")
    print(f"  → {unique_dock} unique compounds exceeded {benchmark_name.title()}'s Avg norm docking")

    print(f"- Compounds outperforming benchmark saved to: {output_path_filtered}")
    print(f"  → {unique_both} unique compounds exceeded {benchmark_name.title()} on both metrics")
