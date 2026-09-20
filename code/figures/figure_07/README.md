# Figure 07

`generate_figure_07.m` reproduces Figure 7, comparing observed water-surface
geometry with theoretical relationships for wavelength, velocity, Froude number,
wave amplitude, and flow depth.

Surface velocity is obtained using KLT-IV.

The figure panels show:

- **(a)** observed wavelength versus surface velocity;
- **(b)** Froude number versus `kh`;
- **(c)** wave amplitude versus `kh`.

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
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_07'))

outputs = generate_figure_07( ...
    'C:\path\to\figure_07_data', ...
    'C:\path\to\figure_07_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the scripts. Generated output is ignored by Git.

The generator writes `Figure7.png` at 600 dpi and `Figure7.pdf` as vector
graphics.

## Software

- MATLAB R2024a or later
