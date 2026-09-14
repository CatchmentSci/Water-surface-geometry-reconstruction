# Figure 08

`generate_figure_08.m` reproduces Figure 8 from the archived real-case
autocorrelation velocity summary.

The panels show:

- **(a)** the deep-water velocity estimate, `U_deep`, against the tracked
  streamwise surface velocity, `U_s`; and
- **(b)** the percentage change in velocity produced by constant, linear and
  power vertical-profile formulations relative to the deep-water estimate,
  plotted against relative depth, `kh = 2*pi*h/lambda`.

Marker colour denotes discharge. Grey bars show the interquartile range in
both coordinates. Cases R12 and R13 are retained in the archived summary for
traceability but are excluded from both panels because their reconstructed
wavelengths were rejected.

## Required archive file

Point the function at the root of the extracted Zenodo archive. It loads:

`videos/outputs/wse_autocorrelation_velocity_method_sensitivity_summary.mat`

Verified SHA-256 checksum:

`301A81689993EB3E315D36B27D654FDC7593C206CB39C840F79B5B511B7B14CA`

The matching human-readable sensitivity table is also deposited as:

`videos/outputs/wse_autocorrelation_velocity_method_sensitivity_summary.csv`

CSV SHA-256:

`C29FF04B094F7BD86C0F9F2AB2BBB770BD17257C59C1AD1D84B9CB4C69E2E7DE`

The MAT-file contains the final 13-case validation and sensitivity tables and
the inclusion mask used for the paper figure. It is a compact post-processing
product: reproducing Figure 8 does not rerun the KLT solver, checkpoint
selection, or wavelength and depth analyses.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_08'))
outputs = generate_figure_08( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_08_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. The generated output is ignored by Git.

The function creates
`wse_autocorrelation_velocity_validation_and_sensitivity.png` at 600 dpi and
the corresponding vector PDF.

## Software

- MATLAB R2024a or later
