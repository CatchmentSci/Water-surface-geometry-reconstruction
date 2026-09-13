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

The analysis and figure-generation workflows use MATLAB. The required MATLAB release and any additional toolboxes will be documented here when the reproducibility package is finalised.

## Repository structure

```text
.
|-- code
|   |-- dependencies    # Supporting functions and third-party dependencies
|   |-- figures         # Scripts used to reproduce paper figures
|   |-- hpc             # High-performance computing scripts and configuration
|   `-- solver          # Inverse-modelling and numerical solver code
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

## Licence

The software in this repository is released under the [MIT License](LICENSE). Data may be subject to additional terms documented alongside the relevant files.
