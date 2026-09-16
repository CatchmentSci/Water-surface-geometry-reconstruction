# Standing-wave wavelength as a function of water depth

`generate_figure_SWW_vs_depth.m` reproduces the figure showing the theoretical
relationship between standing-wave wavelength and water depth for several
surface flow velocities and three assumed vertical velocity profiles.

The figure includes:

- wavelength–depth relationships calculated directly from the theoretical
  wave relationships;
- depth reconstructed from wavelength and surface velocity;
- constant, linear, and power-law vertical velocity profiles; and
- reference lines of constant Froude number.

The calculations do not require external input data.

## Included analysis code

- `generate_figure_SWW_vs_depth.m` is the reproducibility entry point and
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
- Chebfun
- `brewermap` (cmocean/ColorBrewer-compatible colormap utility used for the
  figure colours)

The helper functions `expected_wavelength.m`, `reconstr_depth.m`, and
`expected_kh.m` must be available on the MATLAB path.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'SWW_vs_depth'))

outputs = generate_figure_SWW_vs_depth( ...
    'C:\path\to\figure_output');