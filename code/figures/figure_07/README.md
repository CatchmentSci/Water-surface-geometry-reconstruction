# Figure 07

`generate_figure_7.m` reproduces Figure 7, comparing observed water-surface
geometry with theoretical relationships for wavelength, velocity, Froude number,
wave amplitude, and flow depth.

The figure panels show:

- **(a)** observed wavelength versus surface velocity;
- **(b)** Froude number versus `kh`;
- **(b)** wave amplitude versus `kh`.

## Required data files

Point the function at the root of the dataset. It loads:

`per_transect_initial_accepted.csv`

`Dart_video_statistics.xlsx`

and the profile-summary CSV files contained in:

`real_selected_map_profile_csvs/`

The profile files are matched to the observations using the video date and
time recorded in `Dart_video_statistics.xlsx`.

No additional input files are required.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_07'))

outputs = generate_figure_7( ...
    'C:\path\to\figure_07_data', ...
    'C:\path\to\figure_07_output');