# Figure 03

`generate_figure_03.m` reproduces Figure 3 with alpha `0.85`.

The panels show:

- **(a)** Fr - kh relationship for the constant, linear, and power-function
  profile and for the deep-water approximation;
- **(b)** alpha Fr - kh relationship for the constant, linear, and power-function
  profile and for the deep-water approximation;
  
## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_03'))
outputs = generate_figure_03( ...
    'C:\path\to\figure_03_output');
```

If the second argument is omitted, output is written to an `output` directory
beside the scripts. Generated output is ignored by Git.

The function creates `Figure3.png` at 600 dpi and `Figure3.pdf` as vector
graphics.
  and constructs the first-10-m transect geometry from checkpoint data.

## Software

- MATLAB R2024a or later
- Image Processing Toolbox
