# Water-surface geometry reconstruction for non-contact river flow monitoring through monocular imagery and inverse modeling

This repository contains the data and MATLAB code used to reproduce the analyses and figures presented in the associated research paper by Matthew T. Perks and Giulio Dolcetti.

The study investigates how water-surface geometry and dynamics observed in monocular imagery can be used to estimate river-flow characteristics, including wavelength, amplitude, depth, and velocity.

## Contents

- [About the project](#about-the-project)
- [Prerequisites](#prerequisites)
- [Repository structure](#repository-structure)
- [How to use](#how-to-use)

## About the project

Visible patterns on river surfaces reflect waves, turbulence, and interactions with channel geometry. This project develops methods for reconstructing water-surface geometry from imagery and combining those observations with inverse modelling to support non-contact river-flow monitoring.

## Prerequisites

The operational workflow used MATLAB R2024a on Comet and requires Parallel Computing Toolbox. MATLAB's dependency analyser confirms that the synthetic solver is otherwise self-contained by the files listed below.

## Repository structure

```text
.
|-- code
|   |-- dependencies
|   |   |-- KLT_applyAngleShift.m
|   |   |-- KLT_wrapTo360_centerMedian.m
|   |   |-- LMFnlsq.m
|   |   |-- camera.m
|   |   `-- voxelviewshed.m
|   |-- figures         # Scripts used to reproduce paper figures
|   |-- hpc
|   |   |-- real
|   |   |   |-- klt_analysis_case_lookup.tsv
|   |   |   |-- resubmit_missing_klt_analysis_jobs.sh
|   |   |   |-- run_one_klt_analysis_case.sh
|   |   |   `-- submit_klt_analysis_batch.sh
|   |   `-- syn
|   |       |-- resubmit_missing_klt_jobs.sh
|   |       |-- run_one_klt_case.sh
|   |       |-- run_one_klt_case_ramlog.sh
|   |       `-- submit_klt_syn.sh
|   `-- solver
|       |-- KLT_synthetic_truth_test_v3.m
|       `-- KLT_wse_solver_paths_Apoint_block_jacobi_v5.m
|-- data                # Input and derived data required by the scripts
|-- images              # Images used in this README or other documentation
|-- LICENSE
`-- README.md
```

## How to use

1. Clone or download this repository.
2. Open MATLAB and add the repository to the MATLAB path.
3. Review the relevant script in `code/figures` for its required inputs.
4. Run the script to reproduce the corresponding result or figure.

Detailed instructions, software versions, and the mapping between scripts and paper figures will be added as the repository is populated.

The scripts in `code/hpc` are exact copies of the operational Newcastle University Comet Slurm files and contain environment-specific paths. Update those paths for another system before submission. Synthetic scripts are in `code/hpc/syn`, while real-video scripts are in `code/hpc/real`; their input data and optional checkpoints are supplied through the associated Zenodo dataset.

The current synthetic workflow expects `synthetic_cases.csv` and four Zenodo solver-input files: `solver_inputs_0pt88.mat`, `solver_inputs_1pt10.mat`, `solver_inputs_1pt50.mat`, and `solver_inputs_1pt70.mat`. The RAM-logging worker contains the full 39-case mapping; the standard worker and submission lists are retained exactly as used on Comet for their selected operational runs.

## Licence

The software in this repository is released under the [MIT License](LICENSE). Data may be subject to additional terms documented alongside the relevant files.
