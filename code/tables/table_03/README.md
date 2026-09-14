# Table 03

`generate_table_03.m` reproduces the real-case flow-characteristics table by
combining the deposited hydraulic statistics with seeding densities
recalculated from the 13 solver checkpoints.

## Required archive files

Point the function at the root of the extracted Zenodo archive. It reads:

- `videos/inputs/Dart_video_hydraulic_statistics.csv`;
- `videos/inputs/klt_analysis_case_lookup.tsv`;
- `videos/inputs/sweep_limits.csv`; and
- `videos/outputs/real_wave_wse_pmusic_summary.csv`, which records the
  accepted WSG map for each case; and
- the 13 `<video_identifier>_checkpoint.mat` files in `videos/outputs`.

The hydraulic CSV supplies discharge, flow exceedance, hydraulic depth,
section-averaged velocity, Froude number, and the Reynolds number based on
hydraulic radius. Cases R1-R13 are the records below 124 m3 s-1, ordered by
decreasing discharge.

For each case, seeding density is recalculated as the candidate-point count
within the cleaned, row-restricted adjusted WSG domain divided by that
domain's area. The implementation uses the first valid candidate-count map,
the accepted WSG map reported by the final real-case analysis, a 20-cell
minimum connected area, and the case-specific rows in `sweep_limits.csv`.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'tables', 'table_03'))
outputs = generate_table_03( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\table_03_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The function writes:

- `table_03_flow_characteristics.csv`, including full-precision density
  values and their displayed nearest-integer values; and
- `table_03_flow_characteristics.tex`, containing the complete formatted
  LaTeX table and caption.

## Reproducibility notes

The table uses `ReynoldsNumber_Rh`, not the alternative open-channel value
based on four times the hydraulic radius. The caption's 7-minute interval is
the maximum offset to the nearest supporting flow record; adjacent records
can span 15 minutes.

The deposited hydraulic-statistics CSV is the direct input to this table.
Regenerating that CSV from its original time series additionally requires the
stage/discharge workbook and long-term flow record, which are not included in
the current Zenodo staging archive.

## Software

- MATLAB R2024a or later
- Image Processing Toolbox
