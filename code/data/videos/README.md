# Real-video derived CSV workflow

This folder contains the complete MATLAB calculation chain for rebuilding
the canonical real-video CSV files deposited with the paper. The workflow
starts from the solver checkpoints and velocity-analysis inputs; it does not
use an existing derived CSV as a calculation input.

## Rebuilt products

- `videos/inputs/per_transect_initial_accepted.csv`;
- the 13 `videos/derived/profiles/*_selected_map_profiles_for_real_batch_summary.csv`
  files; and
- `videos/derived/summaries/real_accepted_wsg_summary.csv`.

The profile CSVs contain the selected-map autocorrelation wavelengths,
robust WSG amplitudes, cross-section depths and supporting per-transect
fields. The compact accepted-WSG summary is generated directly from those
13 profile CSVs.

## Required archive inputs

The extracted archive must contain:

- `videos/checkpoints/<identifier>_checkpoint.mat` for R1--R13;
- `videos/inputs/<identifier>_velocity_inputs.mat` for R1--R13;
- `videos/inputs/klt_analysis_case_lookup.tsv`; and
- `videos/inputs/cross_section.csv`.

## Run

In MATLAB R2024a or later:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
archiveRoot = 'C:\path\to\extracted-zenodo-record';
rebuiltRoot = 'C:\path\to\rebuild-check';

addpath(fullfile(repoRoot, 'code', 'data', 'videos'));
outputs = reproduce_video_derived_csvs(archiveRoot, rebuiltRoot);
report = validate_video_derived_csvs(archiveRoot, rebuiltRoot);
```

The rebuild takes several minutes because every checkpoint is reprocessed.
Temporary per-case analysis files are removed automatically after successful
completion. The validator checks the complete table shape, column order,
text fields and numerical values of all 15 canonical CSV files.

## Software

- MATLAB R2024a or later;
- Signal Processing Toolbox (`pmusic`); and
- Statistics and Machine Learning Toolbox (`prctile`).

## File roles

- `reproduce_video_derived_csvs.m` is the public entry point.
- `generate_case_transect_analysis.m` calculates velocity, geometry and
  wave quantities for one checkpoint.
- `export_selected_profile_csvs.m` applies the selector aggregation and
  writes the 13 profile files.
- `export_real_accepted_wsg_summary.m` creates the compact case summary.
- `validate_video_derived_csvs.m` compares rebuilt and deposited products.
- `KLT_select_wave_checkpoint_refined.m` and
  `make_velocity_maps_from_checkpoint.m` are the exact retained calculation
  dependencies used for the publication dataset.

## Relationship to the original analysis scripts

This workflow consolidates the calculations that were previously distributed
across the following analysis scripts:

- `calling_analysis_velocity_export.m`;
- `calling_real_batch_export_selected_profile_csvs.m`;
- `calling_real_batch_analysis_latex_table.m`;
- `batch_first10m_initial_vs_accepted_velocity_scatterplots.m`; and
- `batch_first10m_velocity_validation_and_sensitivity.m`.

Hard-coded machine paths and undocumented manual export steps were replaced
with function arguments and explicit table writers. The selector settings,
quality thresholds, cross-section interpolation, velocity projection,
autocorrelation filtering, robust aggregation and hydraulic equations were
retained. The deposited selector and velocity-map dependencies were copied
without numerical changes; their SHA-256 checksums are respectively
`6e6e506aaf9a7a9df441a6c5e247ded26ee4aaddd4751955bdfbb21fda787569`
and
`6ef451fce8f8ef038ed8676067233854f2f0fd8ffef175f4b7372d48690af637`.

The completed workflow was tested by rebuilding all 15 canonical CSV files
from the deposited checkpoints and inputs. The strict validator found no
differences across 559 table rows.
