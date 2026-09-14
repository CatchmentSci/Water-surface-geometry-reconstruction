# Table B1

`generate_table_b1.m` reproduces Appendix Table B1 from the deposited
synthetic-case autocorrelation summary. The table reports the imposed hydraulic
depth, wavelength, and amplitude together with estimates derived from the
accepted water-surface geometry (WSG) map.

## Required archive file

Point the function at the root of the extracted Zenodo archive. It reads:

- `syn/outputs/synthetic_wave_range_autocorr_summary.csv`.

This is the same accepted-WSG autocorrelation summary used by the Figure 5
workflow. Cases are assigned sequential labels S1-S39 in deposited row order,
rather than relying on legacy labels retained in the source summary.
The archive's `h_m` field is reported as hydraulic depth, `D`, in Table B1.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'tables', 'table_b1'))
outputs = generate_table_b1( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\table_b1_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The function writes:

- `table_b1_synthetic_characteristics.csv`, a clean machine-readable table;
  and
- `table_b1_synthetic_characteristics.tex`, the complete formatted Appendix
  Table B1 and caption.

## Reproducibility notes

The wavelength estimate is the median of streamwise autocorrelation profiles
spaced at 0.5 m in the cross-stream direction. The signed velocity difference
is calculated from the deep-water gravity-wave relation
`U_s = sqrt(g*lambda/(2*pi))`.

Case S19 is the only case whose wavelength and amplitude estimates are
suppressed in the accepted summary; these entries are rendered as `---`.
Large finite differences for S30 and S35 are retained because they are genuine
accepted-summary results.

The LaTeX output requires the `booktabs`, `adjustbox`, and `caption` packages.

## Software

- MATLAB R2024a or later
