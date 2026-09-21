# Synthetic summary reproduction

`reproduce_synthetic_summary.m` rebuilds the accepted autocorrelation summary
used by Figure 5 and Appendix Table B1 directly from the 39 deposited solver
checkpoints. It retains the original map-selection, autocorrelation, robust
amplitude and suppression calculations. PMUSIC is disabled and is not used by
this workflow.

## Run

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
archiveRoot = 'C:\path\to\extracted-zenodo-record';
rebuiltRoot = 'C:\path\to\synthetic-rebuild';

addpath(fullfile(repoRoot, 'code', 'data', 'syn'));
outputs = reproduce_synthetic_summary(archiveRoot, rebuiltRoot);
report = validate_synthetic_summary(archiveRoot, rebuiltRoot);
```

The builder reads the 39 `syn/outputs/*_ckpt.mat` files and writes
`synthetic_wave_range_autocorr_summary.csv` and its MATLAB companion. The
validator checks all fields used by Figure 5 and Table B1, including the
selected map, wavelength, velocity difference, amplitude and suppression flag.

The original analysis was
`calling_syn_batch_analysis_latex_table_simplified.m`. The repository version
replaces hard-coded paths with arguments and removes its legacy PMUSIC and
plot-export branches without changing the accepted autocorrelation results.

## Software

- MATLAB R2024a or later
- Statistics and Machine Learning Toolbox (`prctile`)
