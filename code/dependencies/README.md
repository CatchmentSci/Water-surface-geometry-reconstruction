# Dependencies

| File | Purpose |
|---|---|
| `camera.m` | MATLAB camera-model class used to reconstruct camera objects stored in the solver-input MAT files. |
| `KLT_applyAngleShift.m` | Applies the angular shift used by the synthetic driver. |
| `KLT_wrapTo360_centerMedian.m` | Centres and wraps observed directions before inversion. |
| `LMFnlsq.m` | Nonlinear least-squares routine used by `camera.invproject`; its embedded BSD-style licence and copyright notice are retained. |
| `voxelviewshed.m` | Computes DEM visibility for `camera.prepareDEMInverse`. |

Add this directory to the MATLAB path before running the synthetic-truth workflow. These are byte-for-byte copies of the files used on Comet. MATLAB's dependency analyser identified no additional non-MathWorks source dependencies for the published synthetic workflow.
