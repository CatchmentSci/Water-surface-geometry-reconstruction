# Standing-wave wavelength as a function of water depth

`generate_figure_a2.m` reproduces the figure showing the theoretical
relationship between standing-wave wavelength and water depth for several
surface flow velocities and three assumed vertical velocity profiles.

The figure includes:

- wavelength–depth relationships calculated directly from the theoretical
  wave relationships;
- depth reconstructed from wavelength and surface velocity;
- constant, linear, and power-law vertical velocity profiles; and
- reference lines of constant Froude number.

Profile families are distinguished by colour. Forward wavelength calculations
and inverse depth reconstructions use a common solid line style, while grey
dashed lines show constant-Froude references.
The plot retains the original figure dimensions, axis positions,
limits, logarithmic scaling, and automatic tick placement.
Typography is sized for inclusion at `0.85\textwidth` in the manuscript.

The calculations do not require external input data.

## Included analysis code

- `generate_figure_a2.m` is the reproducibility entry point and
  specifies the output location.
- `plot_SWW_vs_depth.m` performs the theoretical calculations, constructs the
  figure, and exports the results.

The following helper functions provide the underlying calculations:

- `expected_wavelength.m` calculates the standing-wave wavelength from water
  depth and surface velocity by solving the dispersion relationship for the
  wavenumber. It supports constant, linear, and power-law vertical velocity
  profiles.
- `reconstr_depth.m` reconstructs water depth from wavelength and surface
  velocity by solving the corresponding dispersion relationship. Multiple
  mathematical depth solutions may occur; the function returns the first
  solution within the specified depth bounds.
- `expected_kh.m` calculates the dimensionless wavenumber–depth product `kh`
  corresponding to a specified Froude number. The current implementation
  supports the constant-velocity profile.

`expected_wavelength.m` and `reconstr_depth.m` use numerical root finding.
`reconstr_depth.m` uses the Chebfun toolbox to identify the depth solutions
within the prescribed bounds.

## Dependencies

- MATLAB R2024a or later
- Chebfun 5.7.0, matching the version used for this analysis. Download and
  installation instructions are available from the
  [official Chebfun website](https://www.chebfun.org/download/); the archived
  [Chebfun 5.7.0 release](https://github.com/chebfun/chebfun/releases/tag/v5.7.0)
  should be used for exact reproduction. Add the Chebfun root directory to
  the MATLAB path before running the figure generator.

The helper functions `expected_wavelength.m`, `reconstr_depth.m`, and
`expected_kh.m` must be available on the MATLAB path.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_a2'))

outputs = generate_figure_a2( ...
    'C:\path\to\figure_output');
```

The generator writes `FigureA2.png` and `FigureA2.pdf` to the
requested output directory.
