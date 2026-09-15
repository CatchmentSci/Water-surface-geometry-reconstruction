# Water-surface geometry reconstruction for non-contact river flow monitoring through monocular imagery and inverse modeling

This repository contains the MATLAB code used to reproduce the numerical
analyses, figures, and tables supporting the associated research paper by
Matthew T. Perks and Giulio Dolcetti.

The study reconstructs water-surface geometry (WSG) from apparent flow paths
in monocular imagery and uses the recovered surface undulations to estimate
wave properties, surface velocity, and flow depth.

## Repository status

Reproducibility workflows are currently available for:

| Paper item | Repository location | Status |
| --- | --- | --- |
| Figure 1 | `code/figures/figure_01` | Complete; deterministic vector sampling |
| Figure 2 | `code/figures/figure_02` | Final schematic deposited |
| Figure 4 | `code/figures/figure_04` | Complete |
| Figure 5 | `code/figures/figure_05` | Complete |
| Figure 6 | `code/figures/figure_06` | Complete |
| Figure 8 | `code/figures/figure_08` | Complete; initial-planar velocity |
| Figure 9 | `code/figures/figure_09` | Complete; initial-planar velocity |
| Figure 11 (Discussion, depth-bias amplification) | `code/figures/figure_11` | Complete; initial-planar velocity |
| Table 3 | `code/tables/table_03` | Complete |
| Table B1 | `code/tables/table_b1` | Complete |
| Table C1 | `code/tables/table_c1` | Complete |

The remaining figure directories are placeholders for workflows still to be
added. Each completed directory contains its own README with exact archive
inputs, outputs, software requirements, and run instructions.

All velocity-dependent paper figures use
`velocityOutTracked.start.u_streamwise_mps`, calculated by projecting the
tracked paths onto the solver's initial planar water surface. Compact Figure
8 and Figure 9 summaries carrying this choice are version controlled beside
their plotting scripts.

## Data availability

Large inputs and derived products are kept in the associated Zenodo archive,
not duplicated in Git. Point each generator at the root of the extracted
archive, whose top-level data directories are:

```text
syn/
|-- inputs/
`-- outputs/

videos/
|-- inputs/
`-- outputs/
```

The archive contains a SHA-256 manifest and provenance record. Its solver
inputs and checkpoints, together with the repository's operational scripts,
were checked against the read-only Newcastle University Comet working set;
publication-stage summary files document later local post-processing.

Add the final Zenodo DOI and citation here when the record is published.

## Software

The operational solver and HPC workflows used MATLAB R2024a. Individual
figure and table workflows specify their own toolbox requirements; some need
the Image Processing Toolbox or Signal Processing Toolbox. The synthetic HPC
workflow also uses Parallel Computing Toolbox.

## Repository structure

```text
.
|-- code/
|   |-- dependencies/              Shared MATLAB dependencies
|   |-- figures/
|   |   |-- figure_01/ ... figure_11/
|   |   `-- figure_a1/ ... figure_a2/
|   |-- hpc/
|   |   |-- real/                  Real-video Slurm workflow
|   |   `-- syn/                   Synthetic Slurm workflow
|   |-- solver/
|   |   |-- KLT_synthetic_truth_test_v3.m
|   |   `-- KLT_wse_solver_paths_Apoint_block_jacobi_v5.m
|   `-- tables/
|       |-- table_03/
|       |-- table_b1/
|       `-- table_c1/
|-- data/                           Reserved for small repository data
|-- images/                         Documentation images
|-- LICENSE
`-- README.md
```

## Reproducing a figure or table

1. Clone this repository and download/extract the associated Zenodo archive.
2. Open MATLAB and add the selected workflow directory to the path.
3. Call its generator with the archive root and, optionally, an output folder.

For example:

```matlab
repoRoot = 'C:\path\to\Water-surface-geometry-reconstruction';
archiveRoot = 'C:\path\to\extracted-zenodo-archive';

addpath(fullfile(repoRoot, 'code', 'figures', 'figure_05'))
outputs = generate_figure_05(archiveRoot, ...
    fullfile(repoRoot, 'reproduced', 'figure_05'));
```

Generated `output` directories beside figure and table scripts are ignored by
Git. Consult the README in the selected workflow directory before running it.

## Solver and HPC workflows

The shared solver is in `code/solver`, with supporting MATLAB functions in
`code/dependencies`.

The Bash scripts in `code/hpc` are exact copies of the operational Newcastle
University Comet Slurm files. They contain environment-specific paths,
accounts, partitions, log locations, and module names that must be adapted for
another system. Synthetic scripts are in `code/hpc/syn`; real-video scripts
and the matching case lookup are in `code/hpc/real`.

The final synthetic publication set contains 39 cases driven by four archived
solver-input files for hydraulic depths of 0.88, 1.10, 1.50, and 1.70 m. The
RAM-logging worker contains the complete case mapping. The standard worker and
submission lists are retained exactly as used on Comet and record selected
operational runs rather than defining the complete publication set.

The real-video workflow additionally uses the archived case lookup, sweep
limits, solver inputs, observations, and checkpoints. See
`code/hpc/README.md` for the full mapping and portability notes.

## Table conventions

- Table 3 uses the first candidate-count snapshot within the adjusted domain
  accumulated to each accepted WSG map and reports nearest-integer seeding
  density.
- Table B1 uses accepted-WSG autocorrelation estimates for the 39 synthetic
  cases and reports hydraulic depth as `D`.
- Table C1 combines accepted autocorrelation wave estimates with derived
  velocities and initial-planar-velocity deep-branch depth results for stable
  cases R1-R11.

## Licence

The software in this repository is released under the [MIT License](LICENSE).
Data may be subject to additional terms documented with the Zenodo record.
