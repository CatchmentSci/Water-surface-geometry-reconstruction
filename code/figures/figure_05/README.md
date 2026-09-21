# Figure 05

`generate_figure_05.m` reproduces Figure 5 from the archived synthetic-case
summary produced with the autocorrelation wavelength estimator.

The panels compare:

- **(a)** imposed wavelength, `lambda`, with estimated wavelength,
  `lambda_est`; marker colour denotes imposed amplitude; and
- **(b)** imposed amplitude, `A`, with estimated amplitude, `A_est`; marker
  colour denotes imposed wavelength.

The dashed line in each panel is the one-to-one relationship.

## Required archive file

Point the function at the root of the extracted Zenodo archive. It loads:

`syn/outputs/synthetic_wave_range_autocorr_summary.csv`

Verified SHA-256 checksum:

`0DDF8B215EA29B2421E43B069DBFBB016FF78006695DA03E61F1E9E23B47C48C`

The table contains 39 successfully analysed synthetic cases. One case is
flagged for suppression because its stopping criterion did not identify an
amplitude plateau, leaving 38 cases plotted in each panel. Every checkpoint
named by the summary is present in `syn/outputs` in the archive.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_05'))
outputs = generate_figure_05( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_05_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the scripts. The publication-ready PNG and PDF are versioned with the repository.

The function creates `Figure5.png` at 600 dpi and `Figure5.pdf` as vector
graphics.

## Included code

- `generate_figure_05.m` validates the summary data, applies the published
  suppression flag, constructs both panels, and exports the figure.
- `coolwarm2.m` supplies the cool-to-warm diverging colour map used for the
  final paper figure. It is based on the Moreland (2009) colour table.

## Software

- MATLAB R2024a or later
