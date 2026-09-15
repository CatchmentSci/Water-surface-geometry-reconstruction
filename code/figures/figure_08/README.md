# Figure 08

`generate_figure_08.m` reproduces Figure 8 using the streamwise velocity
calculated on the initial planar water surface,
`velocityOutTracked.start.u_streamwise_mps`.

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

## Bundled figure data

The function loads the version-controlled compact summary:

`data/wse_autocorrelation_velocity_method_sensitivity_summary.mat`

The MAT-file contains the final 13-case validation and sensitivity tables,
per-transect initial-planar velocities, and the inclusion mask. Reproducing
Figure 8 does not rerun the KLT solver or wavelength analysis.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_08'))
outputs = generate_figure_08( ...
    '', ... % retained archive argument; bundled data take precedence
    'C:\path\to\figure_08_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. The generated output is ignored by Git.

The function creates
`wse_autocorrelation_velocity_validation_and_sensitivity.png` at 600 dpi and
the corresponding vector PDF.

## Software

- MATLAB R2024a or later
