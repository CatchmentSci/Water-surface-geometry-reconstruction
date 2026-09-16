# Surface velocity as a function of standing-wave wavelength

`generate_figure_Us_vs_SWW.m` reproduces the figure showing the relationship
between surface velocity and standing-wave wavelength for different water
depths and three assumed vertical velocity profiles.

The figure shows:

- the theoretical relationship between surface velocity and wavelength for
  constant, linear, and power-law velocity profiles;
- the deep-water approximation $\sqrt{g/k}$; and
- the corresponding water depths indicated on the secondary y-axis.

No external
input data are required.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'Us_vs_SWW'))

outputs = generate_figure_Us_vs_SWW( ...
    'C:\path\to\figure_output');