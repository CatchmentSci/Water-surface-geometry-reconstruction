# Figure 10

`generate_figure_10.m` reproduces Figure 10, comparing observed water depth
with estimation based on SWA and SWW according to an empirical relationship.

The figure shows observed depth versus estimated depth. The empirical
amplitude--\(kh\) relationship is calibrated from cases with
\(Q>20\ \mathrm{m^3\,s^{-1}}\); the comparison therefore describes
in-sample reconstruction rather than independent validation.

The plotted depth estimate does not depend on velocity. Ancillary velocity
and Froude-number diagnostics returned by the script use KLT-IV surface
velocities.

The plotting style follows Figures 8 and 9: both depth axes span 1--3 m with
equal scaling and matching ticks, uncertainty ranges are drawn beneath the
coloured case markers, the axes use a complete boxed frame, and discharge is
shown using the shared cool-to-warm colour scale at the right of the figure.
Text, marker, and line sizes are pre-scaled for the figure's manuscript
placement at `width=0.6\textwidth` in the 5.5-inch-wide AGU template, so their
displayed sizes match the full-width figures.

## Required data files

These data are supplied by the associated Zenodo repository rather than
duplicated in this code directory. Point the function at the root of that
dataset. It searches subfolders independently for the velocity table and
derived profile-summary directory, and loads:

`per_transect_initial_accepted.csv`

`Dart_video_statistics.xlsx`

and the derived profile-summary CSV files contained in:

`videos/derived/profiles/`

The profile files are matched to the observations using the video date and
time recorded in `Dart_video_statistics.xlsx`.

No additional input files are required.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_10'))

outputs = generate_figure_10( ...
    'C:\path\to\figure_10_data', ...
    'C:\path\to\figure_10_output');
```

The generator writes `Figure10.png` and `Figure10.pdf` to the requested
output directory.
