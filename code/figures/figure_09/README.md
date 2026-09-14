# Figure 09

`generate_figure_09.m` reproduces Figure 9 from the archived real-case
autocorrelation depth-comparison summary.

Panels (a), (b), and (c) compare surveyed cross-section depth, `h`, with the
deep inverse-depth solution, `h_est`, obtained using uniform, linear, and
power velocity-profile formulations, respectively. Marker colour denotes
discharge, grey bars show the interquartile range in both coordinates, and
the dashed line is the one-to-one relationship.

Cases R12 and R13 are retained in the archived table for traceability but are
excluded from the plotted analysis because their reconstructed
autocorrelation wavelengths were rejected.

## Required archive file

Point the function at the root of the extracted Zenodo archive. It loads:

`videos/outputs/wse_autocorrelation_uniform_linear_power_depth_summary.csv`

Verified SHA-256 checksum:

`591FAA058AC6738A184F5D533691C92EC90D5EE544A13B56D42FF242D32BD6BD`

The table contains the final R1-R13 medians, quartiles, asymmetric error
ranges, sample counts, and inclusion flags for both shallow and deep solution
branches. Figure 9 uses the deep branch, while both branches are read to
retain the common axis limits used for the paper's paired depth figures.

This is a compact post-processing product: reproducing Figure 9 does not
rerun the KLT solver, checkpoint selection, wavelength estimation, or
inverse-depth calculations.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_09'))
outputs = generate_figure_09( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_09_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The function creates
`wse_autocorrelation_uniform_linear_power_depth_deep_solutions.png` at 600
dpi and the corresponding vector PDF.

## Software

- MATLAB R2024a or later
