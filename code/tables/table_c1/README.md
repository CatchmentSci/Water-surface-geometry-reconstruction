# Table C1

`generate_table_c1.m` reproduces Appendix Table C1 for the stable real-world
cases R1-R11. It combines accepted autocorrelation wave estimates with
derived velocity and inverse-depth summaries already deposited for Figures 8
and 9.

## Required archive files

Point the function at the root of the extracted Zenodo archive. It reads:

- `videos/outputs/wse_autocorrelation_velocity_method_sensitivity_summary.mat`;
- `videos/outputs/real_wave_wse_pmusic_summary.csv`; and
- `videos/outputs/wse_autocorrelation_uniform_linear_power_depth_summary.csv`.

The MAT-file supplies direct transect-level accepted autocorrelation
wavelengths and case-level velocity medians. The accepted-map summary supplies
the robust WSG amplitude, which does not depend on the wavelength method. The
depth summary supplies the deep-branch median for each velocity profile.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'tables', 'table_c1'))
outputs = generate_table_c1( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\table_c1_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The function writes:

- `table_c1_real_wave_hydraulic_results.csv`; and
- `table_c1_real_wave_hydraulic_results.tex`.

## Reproducibility notes

The generator validates that direct wavelength medians reproduce the
deposited deep-water velocities to the table's three-decimal precision under
`U = sqrt(g*lambda/(2*pi))`, with `g = 9.81 m s^-2`. Tiny full-precision
differences can arise for an even number of transects because the median is
taken before versus after the nonlinear transformation.

R12 and R13 are omitted under the deposited stable-reconstruction inclusion
mask. R5, R7, and R11 have no admissible linear- or power-profile inverse-depth
solution in the deposited Figure 9 summary, so those values are rendered as
`--`. Constant-profile depth values remain available for all R1-R11 cases.

The LaTeX output requires the `booktabs` package.

## Software

- MATLAB R2024a or later
