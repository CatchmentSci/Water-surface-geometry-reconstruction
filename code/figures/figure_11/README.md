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
