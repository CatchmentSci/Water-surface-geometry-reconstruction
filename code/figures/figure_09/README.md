# Figure 09

`generate_figure_09.m` reproduces Figure 9 from depth solutions re-inverted
using KLT-IV surface velocities.

Panels (a), (b), and (c) compare surveyed cross-section depth, `h`, with the
selected inverse-depth solution, `h_est`, obtained using uniform, linear, and
power velocity-profile formulations, respectively. For each transect, the
largest admissible root is selected when multiple roots exist; the sole root
is retained when only one admissible root exists. Case medians can therefore
combine deep-branch and unique-root estimates. Marker colour denotes
discharge, grey bars show the interquartile range in both coordinates, and
the dashed line is the one-to-one relationship.

Cases R12 and R13 are retained in the archived table for traceability but are
excluded from the plotted analysis because their reconstructed
autocorrelation wavelengths were rejected.

## Bundled figure data

The function loads the version-controlled compact summary:

`data/wse_autocorrelation_uniform_linear_power_depth_summary.csv`

The table contains the final R1-R13 medians, quartiles, asymmetric error
ranges, sample counts, and inclusion flags for both shallow and selected
deep-or-sole solutions. Both sets are read to retain the common axis limits
used for the paper's paired depth figures.

`generate_figure_09_depth_summary.m` rebuilds this summary from the
version-controlled Figure 8 transect data and also writes
`data/wse_autocorrelation_depth_root_audit.csv`. The audit table records the
number and value of admissible roots for every case, transect, and velocity
profile. Roots are bracketed on a 0.01 m grid within 0.05--5.00 m; only
subcritical roots are retained. Rebuilding the summary does not rerun the KLT
solver, checkpoint selection, or wavelength estimation.
The rebuild verifies that the Figure 8 data used the shared real-case filter:
finite wavelengths in `0.5 <= lambda <= 7 m` and
`autocorrPeakR >= 0.10` are applied before rejection beyond 3.5 scaled MADs
from the resulting case median (for cases with at least eight remaining
estimates).

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_09'))
depthData = generate_figure_09_depth_summary;
outputs = generate_figure_09( ...
    '', ... % retained archive argument; bundled data take precedence
    'C:\path\to\figure_09_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The plotting function creates `Figure9.png` at 600 dpi and `Figure9.pdf` as
vector graphics.

## Software

- MATLAB R2024a or later
