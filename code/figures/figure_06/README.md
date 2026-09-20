# Figure 06

`generate_figure_06.m` reproduces Figure 6 from the accepted reconstructions
for real cases R1, R3, R5, R7, R9 and R11.

Each panel shows reconstructed water-surface geometry relative to that case's
initial flat surface. All cases are transformed into the R1 streamwise and
cross-stream coordinate frame, retain their manually refined spatial extents,
and use common colour limits. The black boxes mark the regions used for the
paper's detailed comparisons.

## Required archive files

Point the function at the root of the extracted Zenodo archive. It loads the
case-number lookup from:

`videos/inputs/klt_analysis_case_lookup.tsv`

Lookup SHA-256: `A8363B051FC8D739E0DC73E88A92018690AB0147F9A0FBBBC3B7F031CB02C88D`

and the following files from `videos/outputs`:

| Case | Discharge (m3 s-1) | Checkpoint file | SHA-256 |
| --- | ---: | --- | --- |
| R1 | 123 | `devon_dart20181207_12110_checkpoint.mat` | `1BF0575CEA685A211729B1EA05E44616B65FAA2DE060784545FCE39D1FDDC298` |
| R3 | 109 | `devon_dart20181129_13490_checkpoint.mat` | `1841C787FFB77CC747D033A73EA7F1F375D3A3C2E79CAB6EA0E3ED5B27FAA705` |
| R5 | 90 | `devon_dart20181127_14580_checkpoint.mat` | `464F6710ECB7A2CF02EC96119883FC3BFC59E5AA145F1607B1AA0DE0BDDBB7B2` |
| R7 | 70 | `devon_dart20181201_11210_checkpoint.mat` | `C5B7E9798046888E7DDD942054AECA1DCEA6B24DA7D2277C74ED1B585D0293C9` |
| R9 | 51 | `devon_dart20180315_08180_checkpoint.mat` | `5F1DCBE8FFC65861DBF33A45F0C93666E039F7D881B488F0E28E6031467D43A1` |
| R11 | 30 | `devon_dart20180315_17200_checkpoint.mat` | `198703EE83AE07F6D7545698FFC6AE1ACDE4B7AB18226211D0797BE7A8B57C94` |

The case labels are assigned by sorting the lookup table by decreasing
discharge, retaining discharges below 124 m3 s-1, and numbering the remaining
cases sequentially.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_06'))
outputs = generate_figure_06( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_06_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the scripts. Generated output is ignored by Git.

The function creates `Figure6.png` at 600 dpi and `Figure6.pdf` as vector
graphics. The returned structure records the selected checkpoint files,
discharges and accepted map indices.

## Included analysis code

- `generate_figure_06.m` resolves the case numbering, selects each accepted
  map, transforms the grids into the common frame, and exports the figure.
- `KLT_select_wave_checkpoint_refined.m` applies the same accepted-map
  selection used by the real-case batch analysis.

## Software

- MATLAB R2024a or later
- Image Processing Toolbox
