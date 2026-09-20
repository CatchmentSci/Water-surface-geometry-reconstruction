# Surface velocity as a function of standing-wave wavelength

`generate_figure_a1.m` reproduces Appendix Figure A1 showing the relationship
between surface velocity and standing-wave wavelength for different water
depths and three assumed vertical velocity profiles.

The figure shows:

- the theoretical relationship between surface velocity and wavelength for
  constant, linear, and power-law velocity profiles;
- the deep-water approximation $\sqrt{g/k}$; and
- the corresponding water depths indicated on the secondary y-axis.

No external input data are required. The figure is entirely theoretical and
does not read measured velocity data.
The linear and power-profile calculations use \(\alpha=0.85\), wavelengths
span 0.1--100 m, and the five fixed depths span 0.1--10 m. Right-hand depth
labels align with the constant-profile endpoints at the maximum wavelength.
The canvas is calibrated for placement at `width=0.6\textwidth` in
the 5.5-inch AGU manuscript template. At that size, tick and legend text are
approximately 9 pt, axis labels are approximately 10 pt, and the displayed
stroke weights are consistent with the other manuscript figures.

## Run

In MATLAB, from any working directory:

```matlab
repoRoot = 'C:\path\to\repository';
addpath(fullfile(repoRoot, 'code', 'figures', 'figure_a1'))

outputs = generate_figure_a1( ...
    'C:\path\to\figure_output');
```

The generator writes `FigureA1.png` and `FigureA1.pdf` to the requested
output directory.
