# Figure 04

`generate_figure_04.m` reproduces Figure 4 using the archived synthetic-wave
case with imposed amplitude `0.02 m` and wavelength `1.52 m` at the `0.88`
configuration.

The panels show:

- **(a)** angular residual, with the analysed transect band and selected
  central transect overlaid;
- **(b)** reconstructed water-surface elevation relative to the initial flat
  surface; and
- **(c)** imposed and reconstructed detrended water-surface geometry and the
  angular residual along the selected transect.

## Required archive file

Point the function at the root of the extracted Zenodo archive. It loads:

`syn/outputs/syn_solver_inputs_0pt88_case_amp0.02_w1.52_ckpt.mat`

Verified SHA-256 checksum:

`B06FD7D97BAC10EFAF8DAC3327C5828A8B5AD22E42558E5E54C87E87007D17E8`

The checkpoint contains the complete map history, residual maps, grid
coordinates and solver diagnostics required by the plot. No separate solver
input file is needed to reproduce this figure.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_04'))
outputs = generate_figure_04( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_04_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the scripts. Generated output is ignored by Git.

The function creates `Figure4.png` at 600 dpi and `Figure4.pdf` as vector
graphics.

## Included analysis code

- `KLT_plot_refined_wave_solution_three_panel_final.m` constructs the panels
  and exports the figure.
- `KLT_select_wave_checkpoint_refined.m` selects the accepted checkpoint map
  and constructs the first-10-m transect geometry from checkpoint data.

## Software

- MATLAB R2024a or later
- Image Processing Toolbox
