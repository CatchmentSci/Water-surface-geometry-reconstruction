# Figure 10

`generate_figure_10.m` reproduces Figure 10, comparing observed water depth
with estimation based on SWA and SWW according to an empirical relationship.

The figure panels show:

- **(a)** observed depth versus estimated depth.

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

outputs = generate_figure_10( ...
    'C:\path\to\figure_10_data', ...
    'C:\path\to\figure_10_output');