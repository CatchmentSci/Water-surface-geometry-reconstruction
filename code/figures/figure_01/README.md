# Figure 01

`generate_figure_01.m` reproduces Figure 1 using the archived data for the
real-video case `devon_dart20180315_08180` (R9).

## Required archive files

Point the function at the root of the extracted Zenodo archive. It loads:

- `videos/inputs/devon_dart20180315_08180_solver_inputs.mat`
- `videos/checkpoints/devon_dart20180315_08180_checkpoint.mat`

The video itself and the velocity-input MAT file are not required for this
figure because the displayed frame, camera model and flow-vector coordinates
are already stored in the solver-input file.

Verified SHA-256 checksums:

| File | SHA-256 |
| --- | --- |
| solver inputs | `7CA9BCAA76B142B0D499AC51E69CB270FD02F3E2FBF2A3AEEDAD35B55D5857DE` |
| checkpoint | `5F1DCBE8FFC65861DBF33A45F0C93666E039F7D881B488F0E28E6031467D43A1` |

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_01'))
outputs = generate_figure_01( ...
    'C:\path\to\extracted-zenodo-archive', ...
    'C:\path\to\figure_01_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the script. Generated output is ignored by Git.

The function creates:

- `Figure1.png`
- `projected_wse_adjusted_region_with_pixel_angle_vectors_fast.mat`

## Software

- MATLAB R2024a or later
- Image Processing Toolbox
- repository dependencies in `code/dependencies`

The plotting code uses a local Mersenne Twister stream with seed `1` to select
up to 200,000 flow vectors. This makes the selected vector indices identical
on every run without changing MATLAB's global random-number state. The seed and
selected indices are also recorded in the diagnostic MAT file.
