# Figure 11

`generate_figure_11.m` produces the discussion figure showing how a bias in
the tracked surface velocity is amplified by the wavelength-based depth
inversion (Section 4.4).

For a fractional underestimation `epsilon` of the tracked surface velocity,
the constant-profile inversion (Eq. 10a) returns a depth `h_est` satisfying
`tanh(k h_est) = (1 - epsilon)^2 tanh(k h)`, so that

```text
h_est / h = atanh((1 - epsilon)^2 tanh(kh)) / kh
```

The panels show:

- **(a)** the predicted ratio for `epsilon` = 8, 12 and 16 % against
  relative depth `kh = 2*pi*h/lambda`, with the observed ratio for the
  stable R1-R11 cases overlaid. `h` and `h_est` are the paired case medians
  plotted in Figure 9(a) (surveyed cross-section depth and constant-profile
  deep inverse depth over the same transects); `lambda` is the accepted
  autocorrelation wavelength of Table C1. Grey bars show the interquartile
  ranges of `h` (mapped to `kh`) and of `h_est` (divided by the median `h`).
  Marker colour denotes discharge.
- **(b)** the same ratio as a function of `epsilon` for `kh` = 1, 3 and 8,
  bracketing the field range, with the 8 % reference bias (the median
  `U_deep`/`U_s` offset of Figure 8a) marked.

Cases R12 and R13 are retained in the archived summary but excluded from
panel (a) under the shared stable-reconstruction mask.

## Required archive files

Point the function at the root of the extracted Zenodo archive. It loads:

`videos/outputs/wse_autocorrelation_uniform_linear_power_depth_summary.csv`

SHA-256:

`591FAA058AC6738A184F5D533691C92EC90D5EE544A13B56D42FF242D32BD6BD`

`videos/outputs/table_c1_real_wave_hydraulic_results.csv`

SHA-256:

`95A58E56DFF92C685E44E55380762AF2D4A246FE0262C64F0EEE493D05D83457`

The function checks that `depthConstant_m` in Table C1 matches the deep
uniform-profile median in the Figure 9 summary and stops if the two archive
products disagree.

This is a compact post-processing product: reproducing the figure does not
rerun the KLT solver, checkpoint selection, wavelength estimation, or
inverse-depth calculations.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_11'))
outputs = generate_figure_11( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_11_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

Optional name-value arguments override the curves drawn:

```matlab
outputs = generate_figure_11(archiveRoot, outputFolder, ...
    "EpsilonValues", [0.08 0.12 0.16], ...   % panel (a) curves
    "KhValues", [1 3 8], ...                 % panel (b) curves
    "ReferenceEpsilon", 0.08);               % panel (b) marker
```

The function creates
`wse_autocorrelation_depth_bias_amplification.png` at 600 dpi, the
corresponding vector PDF, and a CSV of the plotted values with the predicted
ratio at each case and the velocity bias that would reproduce each observed
ratio exactly (`impliedVelocityBias`).

## Software

- MATLAB R2024a or later

## Velocity-source variant

`generate_figure_11_velocity_variant.m` rebuilds the observed points of
panel (a) from a transect-level re-inversion driven by a chosen tracked
surface velocity, so that the figure can be produced with the velocity
projected onto the accepted water-surface map instead of the velocity
deposited with the archived depth summary. Panel (b) is unchanged.

It needs, in addition to the archive files above,

- `videos/outputs/wse_autocorrelation_velocity_method_sensitivity_summary.mat`
  (Figure 8 per-transect tables: accepted autocorrelation wavelength and the
  `commonValid` mask), and
- `checkpoint_batch_initial_vs_accepted_velocity_summary.mat` written by
  `batch_first10m_initial_vs_accepted_velocity_scatterplots.m`
  (`initialVelocityTrackedMedian_mps` and `acceptedVelocityTrackedMedian_mps`
  per transect).

On the originating workstation, the second file is available directly at:

```text
D:\OneDrive - Newcastle University\Documents - WSE Project\General\Dart\Videos\Inputs\batch_first10m_initial_vs_accepted_velocity_outputs\checkpoint_batch_initial_vs_accepted_velocity_summary.mat
```

This path is the optional default in the function. Pass an explicit
`batchMatFile` on another computer or after moving the file.

Per transect the constant-profile inversion (Eq. 10a) has the single
admissible root `h_est = atanh(U^2 k / g) / k`, defined only where
`U^2 k / g < 1`; case values are paired medians of surveyed depth and
`h_est` over the `commonValid` transects with an admissible root, as for
Figure 9(a).

```matlab
outputs = generate_figure_11_velocity_variant(archiveRoot, ...
    'C:\path\to\checkpoint_batch_initial_vs_accepted_velocity_summary.mat', ...
    outputFolder, "VelocitySource", "accepted");   % | "initial" | "archived"
```

On the originating workstation, the default direct MAT-file link can be used
by supplying only the archive root (and accepting the default output folder):

```matlab
outputs = generate_figure_11_velocity_variant(archiveRoot);
```

`"archived"` reproduces the deposited Figure 9/11 values without
re-inversion and is the reference the other two sources are compared to.
The function prints, per case, the maximum absolute difference between the
archived per-transect velocity and each batch velocity, so the provenance
of the deposited figure is checked rather than assumed. `"MaxDepth_m"`
discards transect inversions deeper than a threshold (default `Inf`).

Outputs: `wse_autocorrelation_depth_bias_amplification_<source>_velocity.{png,pdf}`,
`figure_11_velocity_source_comparison_by_case.csv` (kh, ratio, IQR bars and
implied velocity bias for all three sources, plus the ratio differences to
the archived values) and `figure_11_velocity_source_comparison_by_transect.csv`.
